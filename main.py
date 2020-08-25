from Scripts import (dataset, experiment, settings)
from pathlib import Path
import json

kmers = settings.KMERS
kraken_filter = settings.KRFILTER

def main():
    Path("Datasets").mkdir(exist_ok=True)
    Path("Database").mkdir(exist_ok=True)

    data = json.load(open("Settings/datasets.json"))
    deamination = json.load(open("Settings/deamination.json"))

    datasets = []
    regimes = list(deamination.keys())

    #Create the datasets
    for record in data:
        datasets.append(record)

        dataset.make(record, data[record], deamination)

    #Now run the experiments
    for ds in datasets:
        for regime in set(regimes):
            for kmer in set(kmers):
                for kf in kraken_filter:
                    experiment.main(ds, regime, kmer, kf)


if __name__ == "__main__":
    main()
