from Scripts import dataset
import os
from pathlib import Path
import json


def main():
    datadir = Path("Datasets")
    if not datadir.exists():
        datadir.mkdir()
    
    if not Path("Database").exists():
        Path("Database").mkdir()

    data = json.load(open("Settings/datasets.json"))
    deamination = json.load(open("Settings/deamination.json"))
    
    for record in data:
        dataset.make(record, data[record], deamination)


if __name__ == "__main__":
    main()
