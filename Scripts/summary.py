from pathlib import Path
from . import settings
import pysam
import requests
from bs4 import BeautifulSoup

taxonomy = {}


class Sequence:
    def __init__(self, description, dataset, regime, kmer, kraken_filter):
        self.description = description
        self.dataset = dataset
        self.regime = regime
        self.kmer = kmer
        self.kraken_filter = kraken_filter
        self.accession = None
        self.original_family = None
        self.assigned_family = None
        self.extracted = False
        self.aligned = False
        self.bedfiltered = False

    def get_accession_from_description(self):
        return self.description.split("|")[0]

    def get_taxonomy_from_accession(self):
        if self.accession is None:
            self.accession = self.get_accession_from_description()
        if self.accession in taxonomy:
            self.original_family = taxonomy[self.accession]
            return True
        else:
            soup = BeautifulSoup(
                requests.get(
                    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=taxonomy&id=' +
                    self.accession
                ).text, features="html.parser"
            )
            try:
                tax_id = soup.find("linksetdb").find("id").contents[0]
                soup2 = BeautifulSoup(
                    requests.get(
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tax_id}"
                    ).text, features="html.parser"
                )
                for taxon in soup2.find("lineageex").find_all("taxon"):
                    if taxon.rank.string == "family":
                        taxonomy[self.accession] = taxon.scientificname.string
                        self.original_family = taxon.scientificname.string
                        return True
                    else:
                        self.original_family = "Unknown"
                        return False
            except AttributeError:
                self.original_family = "Unknown"
                return False


def extract_records_from_bam(samfile):
    try:
        return [x.split("\t")[0] for x in open(samfile, "r")]
    except:
        pysam.index(str(samfile))
        return [x.query_name for x in pysam.AlignmentFile(samfile, check_sq=False)]


def summarize_megan(location, records, ds_id, regime):
    print(f"Summarize MEGAN: Dataset{ds_id}, Regime{regime}")
    rec_seq_dict = {rec: Sequence(rec, ds_id, regime, "MEGAN", 0) for rec in records}
    run_dir = location.joinpath(f"MEGAN")

    for family in run_dir.joinpath("out").joinpath("blast").iterdir():
        extracted_bam = [x for x in family.glob("*.bam")]
        if len(extracted_bam) < 1:
            continue
        for record in extract_records_from_bam(extracted_bam[0]):
            rec_seq_dict[record].extracted = True
            rec_seq_dict[record].assigned_family = family.name
        for aligned_bam in [x for x in family.joinpath("aligned").glob("*.bam")]:
            for record in extract_records_from_bam(aligned_bam):
                rec_seq_dict[record].aligned = True

    return list(rec_seq_dict.values())


def summarize_kraken(location, records, ds_id, regime, kmer, krf):
    print(f"Summarize Kraken: Dataset{ds_id}, Regime{regime}, Kmer{kmer}, Krakenfilter {krf}")
    print(f"Records in kraken: {len(records)}")
    rec_seq_dict = {rec: Sequence(rec, ds_id, regime, kmer, krf) for rec in records}
    run_dir = location.joinpath(f"Kmer{kmer}").joinpath(f"kraken_filter_{krf}")
    print(f"ref_seq_dict: {len(rec_seq_dict)}")
    for family in run_dir.joinpath("out").iterdir():
        extracted_bam = [x for x in family.glob("*.bam")]
        if len(extracted_bam) < 1:
            continue
        for record in extract_records_from_bam(extracted_bam[0]):
            rec_seq_dict[record].extracted = True
            rec_seq_dict[record].assigned_family = family.name
        for aligned_bam in [x for x in family.joinpath("aligned").glob("*[!_deduped].bam")]:
            for record in extract_records_from_bam(aligned_bam):
                rec_seq_dict[record].aligned = True
        for bedfiltered_bam in [x for x in family.joinpath("bedfiltered").glob("*.bam")]:
            for record in extract_records_from_bam(bedfiltered_bam):
                rec_seq_dict[record].bedfiltered = True

    return list(rec_seq_dict.values())


def main(ds_id, regime, dataset_dir="Datasets", savedir="Summary", experiment_dir="Experiments"):
    """
    for the given dataset_id and the given regime identifier, summarize the results of this run and save it in the
    summary-directory
    Args:
        ds_id:              the number of the dataset
        regime:             the mode
        dataset_dir:        the pointer to the original dataset (default: Dataset)
        savedir:            the directory to save the summary in (default: Summary)
        experiment_dir:     the pointer to the Experiments (default: Experiments)

    Returns:
        True on success, writes files to the system
    """
    savedir = Path(savedir)
    dataset_dir = Path(dataset_dir)
    experiment_dir = Path(experiment_dir)

    kmers = settings.KMERS
    kfilter = settings.KRFILTER

    all_seqs = []

    records = [x.split("\t")[0] for x in open(dataset_dir.joinpath(f"Dataset_{ds_id}_{regime}.sam"), "r")]
    print(f"Records:{len(records)}")
    location = experiment_dir.joinpath(f"Dataset{ds_id}").joinpath(f"Regime{regime}")

    for kmer in kmers:
        if kmer == "MEGAN":
            seqs = summarize_megan(location, records, ds_id, regime)
            print(len(seqs))
            all_seqs.extend(seqs)
        else:
            for krf in kfilter:
                seqs = summarize_kraken(location, records, ds_id, regime, kmer, krf)
                print(len(seqs))
                all_seqs.extend(seqs)

    savedir.mkdir(exist_ok=True, parents=True)
    with open(savedir.joinpath(f"Summary_Dataset{ds_id}_Regime{regime}.tsv"), "w") as outfile:
        outfile.write(
            "Description\tAccession\tDataset\tRegime\tKmer\tKrakenfilter\tFamily\t"
            "Assigned_Fam\tExtracted\tAligned\tBedfiltered\n"
        )
        for seq in all_seqs:
            outfile.write(
                f"{seq.description}\t{seq.get_accession_from_description()}\t{seq.dataset}\t{seq.regime}\t{seq.kmer}\t{seq.kraken_filter}\t"
                f"{seq.original_family}\t{seq.assigned_family}\t{seq.extracted}\t{seq.aligned}\t{seq.bedfiltered}\n"
            )
    return True
