from pathlib import Path
from Bio import (Entrez, SeqIO)
from Bio.Seq import Seq
import requests
import gzip
import io
import sys
from bs4 import BeautifulSoup
from Scripts import (chunk_genome, settings)
import random
import itertools
import os
import json


def create_samfiles(directory):
    """
    Creates a samfile from the fastafiles within a directory.
    Args:
        directory:  The path to the directory with the fasta-files
    Returns: True on success
    """
    for infile in directory.glob("*.fas"):
        with open(directory.joinpath(f"{infile.stem}.sam"), "w") as outsam:
            fasta = SeqIO.parse(open(infile), "fasta")
            for record in fasta:
                quality = "]" * len(record.seq)
                line = f"{record.description}\t4\t*\t0\t0\t*\t*\t0\t0\t" \
                       f"{record.seq}\t{quality}\tXI:Z:TAATCAT\tYI:Z:]]]]]]]\tXJ:Z:TCATGGT\t" \
                       f"YJ:Z:]]]]]]]\tRG:Z:Simdata\tNM:i:0\n"
                outsam.write(line)
        infile.unlink()
    with open(directory.joinpath("indexfile.tsv"), "w") as indexfile:
        indexfile.write("#Index\tlibID\tprimer_p7\tprimer_p5\nSimdata\t368\t34\trandom\n")
    return True


def merge_snippets(directory, ds_id, regimes):
    """
    This process groups all genome-snippets from mammals and contaminants by
    Dataset and Regime identifier and merges them into one Dataset (each)
    Args:
        regimes: a list of regimes
        ds_id: the identifier of the dataset
        directory: The directory in which the chunked and deaminated genomes are stored

    Returns True after the process is done
    """

    for reg in list(regimes):
        outpath = directory / Path(f"Dataset_{ds_id}_{reg}.fas")
        with open(outpath, "w") as outfile:
            for path in Path(directory).glob(f"{ds_id}_{reg}*"):
                for record in SeqIO.parse(open(path), "fasta"):
                    SeqIO.write(record, outfile, "fasta")
                Path(path).unlink(missing_ok=True)
    return True


def deaminate_genome(inpath, outpath, positions=3, prop=0.5, inner=False, prop2=0):
    """
    Deaminate a Sequence
    Args:
        inpath:       The chunked genome [path]
        outpath:    The name of the output-file
        positions:   Number of positions checked for deamination [int]
        prop:        Chance for a C to T substitution [float]
        inner:       Also deaminate inner Cs? [bool]
        prop2:      With that chance [float]

    Returns the deaminated Sequence
    """
    infile = SeqIO.parse(open(inpath), "fasta")
    with open(outpath, "w") as outfile:
        for record in infile:

            new_seq = [x for x in str(record.seq)]

            for n in range(positions):
                new_seq[n] = "T" if new_seq[n] == "C" and random.random() < prop else new_seq[n]
                new_seq[-n - 1] = "T" if new_seq[-n - 1] == "C" and random.random() < prop else new_seq[-n - 1]

            if inner:
                for n in range(positions + 1, len(new_seq) - positions):
                    new_seq[n] = "T" if new_seq[n] == "C" and random.random() < prop2 else new_seq[n]
            record._set_seq(Seq("".join(new_seq)))
            SeqIO.write(record, outfile, "fasta")
    return True


def download_genome(acc, savedir="Database"):
    """Download a genome based on the NCBI Accession number
    Args:
        acc:	    NCBI Accession Number of the genome you want
        savedir:    Optional: specify the location/directory for the file

    Returns True on Success
    """
    soup = BeautifulSoup(requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=taxonomy&id=' + acc).text,
                         features="html.parser")
    try:
        seq_id = soup.find("idlist").find("id").contents[0]
    except AttributeError:
        print(f"Error: Accession Number {acc} not found", file=sys.stderr)
        return False

    Entrez.email = settings.EMAIL
    infile = Entrez.efetch(db="Nucleotide", id=seq_id, rettype="fasta")
    infile2 = io.BytesIO(bytes(infile.read(), encoding="utf-8"))

    Path(savedir).mkdir(exist_ok=True, parents=True)
    savedir = str(Path(savedir) / f"Sequence_{acc}.fas.gz")

    handle = gzip.open(savedir, "wb")
    handle.write(infile2.read())
    handle.close()
    return True


def download_contaminants(savedir="Database"):
    outfile = gzip.open(Path(savedir).joinpath("Sequence_Contamination.fas.gz"), "wb")
    acclist = [x.replace("\n", "") for x in open("Settings/contamination_list.txt")]
    Entrez.email = settings.EMAIL
    chromosomes = len(acclist)
    for i in range(len(acclist)):
        bacterium = acclist[i]
        try:
            handle = Entrez.esearch(db="nucleotide", term=bacterium)
            record = Entrez.read(handle)
            handle.close()
            acc = record["IdList"][0]
        except IndexError:
            chromosomes -= 1
            continue
        print(f"Create contamination Dataset. Download Sequence {i+1}/{len(acclist)}", end="\r", file=sys.stderr)
        if download_genome(acc, savedir):
            bact = Path(savedir).joinpath(f"Sequence_{acc}.fas.gz")
            # test, if the bacterial genome is of a good size
            size_test = SeqIO.parse(gzip.open(bact, "rt"), "fasta")
            for record in size_test:
                if len(record) < 120:
                    bact.unlink()
                    chromosomes -= 1
                    print(f"Record {record.id} removed, genome too small", end="\n", file=sys.stderr)
            if bact.exists():
                infile = gzip.open(bact, "rb")
                outfile.write(infile.read())
                infile.close()
                bact.unlink()
    outfile.close()
    return chromosomes, True


def make(ds_id, data, deamination, outdir="Datasets"):
    """Create a Dataset based on the
    json file provided.

    Args:
        deamination: Dict with the Deamination Settings
        ds_id: 	Identifier of the Dataset
        data: 	Dict with the following structure:
                {Accession:{Reads:1000,deamination:{B:0,C:0}}}, Accession: {...}}
        outdir: Directory where to save the created Datasets (default: Datasets)

    Writes all the files to the file-system and returns True on success
    """
    database = Path("Database")
    datasets = Path(outdir)
    datasets.mkdir(exist_ok=True)
    bact_chromosomes = 70

    for acc in data:
        if f"Sequence_{acc}.fas.gz" not in [x.name for x in database.glob("*.gz")]:
            if acc != "Contamination":
                while True:
                    if download_genome(acc, savedir=database):
                        break
            else:
                while True:
                    check, bact_chromosomes = download_contaminants()
                    if check:
                        break
        print(f"Processing File Sequence_{acc}.fas.gz", end="\r", file=sys.stderr)
        infile = database.joinpath(f"Sequence_{acc}.fas.gz")
        outfile = datasets.joinpath(f"Sequence_{acc}_chunked.fas")
        num_seq = data[acc]["Reads"]
        length = Path("Settings/read_length_dist.tsv")
        chromosomes = 1 if acc != "Contamination" else bact_chromosomes
        minlen = 35
        maxlen = 100
        unif = 0.01

        check = chunk_genome.main(infile, outfile, num_seq, length, chromosomes, minlen, maxlen, unif)
        if check:
            for setup in data[acc]["Deamination"]:
                save_file = datasets.joinpath(f"{ds_id}_{setup}_{acc}_chunked.fas")

                deaminate_genome(outfile, save_file, positions=deamination[setup]["Ends"]["Positions"],
                                 prop=deamination[setup]["Ends"]['Probability'],
                                 inner=True, prop2=deamination[setup]["Center"]['Probability'])
        outfile.unlink()

    merge_snippets(datasets, ds_id, data[acc]["Deamination"].keys())
    create_samfiles(datasets)
    return True
