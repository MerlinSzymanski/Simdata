from pathlib import Path
from Bio import Entrez

from . import chunk_genome
import re, os, sys, shutil, random
import json, requests, gzip
from bs4 import BeautifulSoup


def download_genome(acc):
	"""This function downloads a genome from NCBI based on the Accession number
	Args:
		acc:	NCBI Accession Number of the genome you want

	Returns:
				True on success
	"""
	soup = BeautifulSoup(requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=taxonomy&id='+acc).text)
	seq_id = soup.find("idlist").find("id").contents[0]

	Entrez.email = "merlin.szymanski@gmail.com"
	infile = Entrez.efetch(db="Nucleotide", id=seq_id, rettype="fasta")

	#TODO: GZIP CRASHES
	with gzip.open(f"Database/Sequence_{acc}.fas.gz", "wb") as out:
		data = bytes(infile.readlines(), encode="UTF-8")
		out.write(data)

	return True


def make(id, data):
	""" The main function of this module. It creates a Dataset based on the 
	json file provied.
	
	Args: 	
		id: 	Identifier of the Dataset [int]
		data: 	Dict with the following structure:
				{Accession:{Reads:1000,Deamination:{B:0,C:0}}}, Accession: {...}}
	Returns: 
		True on success
	"""

	dataset = Path("Datasets"/f"Dataset{id}")
	#Todo create 3 temp dirs A,B and C

	for accession in data:
		#1. get the right mt genome
		#2. chunk it

		#for regime in deamination
			#for read in chunked_genome
				#3. deaminate and put into corresponding directory

	#bundle everything together
		return True