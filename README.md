# Simdata #

This Script:
 - Creates simulated Datasets of ancient DNA
 - Runs them on the Pipeline(s) 
 - WIP: Summarizes the Results  

## Getting started ##

Make sure, everything one needs is installed/prepared:
- Nextflow pipeline on hand
- Snakmake pipeline on hand
- Kraken installed  
- Kraken Databases indexed
- Genomes for BWA prepared

If you are at the **MPI Eva**, you can ignore the things above and \
reproduce my work by just running:

```
git clone https://www.github.com/MerlinSzymanski/Simdata
cd Simdata
conda env create -f environment.yml
conda activate simdata
python3 main.py
```

If you are somewhere else... well, shit

## Change Things ##
#### 1. The Datasets ####
**The endogeneous DNA**

The Datasets are defined in the json file in ```Settings/datasets.json```
what you need to specify:  
- a Dataset number
- the accession number from NCBI
- the number of reads
- for the modes specified in ```Settings/deamination.json```, if they apply

I kow, this is hard to read. Maybe I'll write a GUI for that in the future

The raw downloaded and gziped genomes are saved in the ```Database``` directory, 
the created datasets are stored in the ```Datasets``` folder


**The Deamination**

For the modes specified in ```Settings/deamination.json``` one can specify:
- Mode identifier (usually A,B,C)
- The number of positions on the ends of the sequence
- The chance for a C-T substitution
- If inner C's are also substituted
- The chance for that to happen

**The Contamination**

```Settings/contamination_list.txt```\
The contamination is just a list of 110 genomes (Bacteria, Archae). The genomes
are searched by name and downloaded from NCBI.

**Read length distribution**

```read_length_dist.tsv```\
This file specifies the distribution of read lengths for the genomes when the 
genomes are chunked into pieces 

#### 2. The Experiment ####

Most of the experimental settings one can find in ```Scripts/settings.py```

Please specify there the location of your pipelines, a list of kmers to test, 
the location of the genomes etc.
