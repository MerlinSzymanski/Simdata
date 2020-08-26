FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml .

RUN conda env create -f environment.yml

SHELL ["conda", "run", "-n", "simdata", "/bin/bash", "-c"]

COPY . .

ENTRYPOINT ["conda", "run", "-n", "simdata", "python3", "main.py"]

