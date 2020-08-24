## Simdata ##

### Create the Datasets ###

this project allows a completely reproducible workflow
To create the Datasets used in the analysis run:

* conda env create -f environment.yml
* conda activate simdata
* python3 main.py

The raw gziped genomes are saved in the "Database" folder, the 
datasets are stored in the "Datasets" directory

To change the **content** of the datasets use the files in the "Settings" directory
* datasets.json: here, for each accession numbers for the "endogenous" DNA, one can set up the number of reads and which modifications apply
* deamination.json: Here you can modify the modifications
* contamination\_list.json: A list of bacterial species used to create the contamination 
* read\_length\_dist: The distribution of read lengths from chargirskaya cave used for sampling

please go into Scripts/settings.py and change the EMAIL used for NCBI Entrez

