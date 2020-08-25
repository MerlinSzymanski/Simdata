from pathlib import Path
import subprocess
import os
from . import settings


def run_megan(dataset, index, wdir):
    cwd = os.getcwd()
    dataset = dataset.absolute()
    index = index.absolute()
    wdir.mkdir(exist_ok=True, parents=True)
    os.chdir(wdir.absolute())

    orders1 = ["snakemake", "--jobs", "48", "-s", settings.SNAKEFILE1, "--config",
               f"bamfile={str(dataset)}", f"byfile={str(index)}"]
    orders2 = ["snakemake", "--jobs", "48", "-s", settings.SNAKEFILE2]
    orders3 = ["snakemake", "--jobs", "48", "-s", settings.SNAKEFILE3]

    print("Run: ", " ".join([str(x) for x in orders1]))
    subprocess.run(orders1)
    print("Run: ", " ".join([str(x) for x in orders2]))
    subprocess.run(orders2)
    print("Run: ", " ".join([str(x) for x in orders3]))
    subprocess.run(orders3)

    os.chdir(cwd)
    return True


def run_kraken(dataset, index, wdir, kmer, kraken_filter):
    dataset = dataset.absolute()
    index = index.absolute()

    nextflow = settings.NEXTFLOW
    kmerdir = settings.KMERDIR
    genomes = settings.GENOMES

    cwd = os.getcwd()
    wdir.mkdir(exist_ok=True, parents=True)
    os.chdir(wdir.absolute())

    orders = ["nextflow", "run", "-resume", nextflow, "--rg", str(index), "--bam",
              str(dataset), "--db", f"{kmerdir}/Mito_RefSeqRel97_Kmer{kmer}", "--genome",
              genomes]

    if kraken_filter > 0:
        orders.extend(["--krakenfilter", kraken_filter])

    print(f"Running:{' '.join([str(x) for x in orders])}")
    subprocess.run(orders)

    os.chdir(cwd)
    return True


def main(ds_id, regime, kmer, kraken_filter):
    dataset = Path(f"Datasets/Dataset_{ds_id}_{regime}.sam")
    index = Path(f"Datasets/indexfile.tsv")

    if kmer != "MEGAN":
        wdir = Path(f"Experiments/Dataset{ds_id}/Regime{regime}/Kmer{kmer}/kraken_filter_{kraken_filter}")
        run_kraken(dataset, index, wdir, kmer, kraken_filter)

    else:
        wdir = Path(f"Experiments/Dataset{ds_id}/Regime{regime}/MEGAN")
        run_megan(dataset, index, wdir)

    return True
