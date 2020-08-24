from pathlib import Path
import os
import settings

def run_megan(dataset, wdir, index):
    os.mkdir("MEGAN")
    orders1 = ["snakemake","--jobs", "48", "-s", "/mnt/scratch/merlin/MEGAN_pipeline/metagen.p1.SnakeFile",
                    "--config", f"bamfile={dataset.absolute()}",  f"byfile={index.absolute()}"]
    orders2 = ["snakemake","--jobs", "48","-s", "/mnt/scratch/merlin/MEGAN_pipeline/metagen.p2.SnakeFile"]
    orders3 = ["snakemake","--jobs", "48","-s", "/mnt/scratch/merlin/MEGAN_pipeline/metagen.p3.SnakeFile"]

    print("Run: ", " ".join(orders1))
    subprocess.run(orders1)
    print("Run: ", " ".join(orders2))
    subprocess.run(orders2)
    print("Run: ", " ".join(orders3))
    subprocess.run(orders3)


def run_kraken(dataset, wdir, kmer, index):
    nextflow = settings.NEXTFLOW
    kmerdir = settings.KMERDIR
    genomes = settings.GENOMES
    
    cwd = os.getcwd() 
    os.chdir(wdir.joinpath(f"Kmer{kmer}").mkdir(exists_ok=True))
    
    orders = ["nextflow", "run", "-resume", nextflow, "--rg", index.absolute(), "--bam", 
            dataset.absolute(), "--db", f"{kmerdir}/Mito_db_Kmer{kmer}", "--genome", 
            genomes]
    print(f"Running:{' '.join(orders)}")
    subprocess.run(orders)

    os.chdir(cwd)
    return True

def main(ds_id, regime):
    dataset = Path(f"Datasets/Dataset_{ds_id}_{regime}")
    index = Path(f"Datasets/indexfile.tsv")
    wdir = Path(f"Experiments/Dataset{ds_id}/Regime{regime}/").mkdir(exists_ok=True, parents=True)
    for kmer in settings.KMERS:
        run_kraken(dataset, wdir, kmer, index)
    
    run_megan()
    return True

