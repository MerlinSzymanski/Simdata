import gzip
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import random
import sys
import re
import argparse
import pysam

random.seed()


def sort_recs(recs, split_char="|"):
    def sort_func(item):
        # we sort according to chromosomes
        # our header is >16|69694935_57|Nea ...
        # use split to access the first element
        return item.id.split(split_char)[0]
    return recs


def mutate_unif(sequence, unif):
    """
    Create uniform distributed random mutations in the sequence
    Args:
        sequence:       The DNA-Sequence to be mutated
        unif:    Number between 0 and 1, chance for a mutation to occur

    Returns the mutated Sequence
    """
    same_nuc = set("ACGT")

    choices = np.random.random(len(sequence)) < unif
    new_seq = [''] * len(sequence)

    for idx, nc in enumerate(sequence):
        mutate = choices[idx]
        if mutate:
            new_seq[idx] = random.choice(tuple(same_nuc.difference(nc)))
        else:
            new_seq[idx] = nc
    return Seq(''.join(new_seq))


def chunk_fast(record, n_samples, unif=None, len_distrib=None, minlength=35, maxlength=100):
    try:
        positions = random.sample(range(0, len(record) - maxlength), n_samples)
    except:
        # sample too small sample with replacement
        positions = [random.choice(range(0, len(record) - maxlength)) for _ in range(n_samples)]

    if len_distrib:  # we gave a file with distribution per length
        p = read_len_distrib(len_distrib, minlength, maxlength)
        length = np.random.choice(np.arange(minlength, maxlength + 1), n_samples, p=p)
    else:
        length = [random.choice(range(minlength, maxlength + 1)) for k in range(n_samples)]

    all_samples = []
    for pos, size in zip(positions, length):
        while 'N' in record[pos:pos + size]:
            pos = random.choice(range(0, len(record) - maxlength))
            if len_distrib:
                size = np.random.choice(
                    np.arange(minlength, maxlength + 1), 1, p=p)[0]
            else:
                size = random.choice(range(minlength, maxlength + 1))
        sample = record[pos:pos + l]

        if unif:
            sample.seq = mutate_unif(sample.seq.tomutable(), unif)
        all_samples += [(sample, pos)]
    return all_samples


def read_len_distrib(filename, minlen=35, maxlen=100):
    from csv import reader
    with open(filename, 'r', newline='') as csvfile:
        csvreader = reader(csvfile, delimiter='\t')
        first_row = csvreader.__next__()
        read_distrib = [0] * int(first_row[1]) + [int(first_row[0])]
        read_distrib.extend([int(row[0]) for row in csvreader])
        # print(read_distrib)
    if len(read_distrib) < maxlen:
        print("Read length distribution provides values only until",
              len(read_distrib), "; expected:", maxlen)
    s = sum(read_distrib[minlen:maxlen + 1])
    return [r / s for r in read_distrib[minlen:maxlen + 1]]


def estimate_read_distribution(file_in, num_seq, n_chromosomes=None):
    """
    This method takes the total number of reads and the number of chromosomes and gives back
    the number of samples per chromosome

    Returns: A list of number of reads
    """
    try:
        with pysam.FastaFile(file_in) as fa:
            sample_n = n_chromosomes if n_chromosomes else fa.lengths

            full_size = sum(fa.lengths[:sample_n])
            size_percent = [int(s / full_size * num_seq) + 10 for s in fa.lengths[:sample_n]]
            return size_percent

    except OSError:  # fasta is not bgzip'd do a naive fallbac
        try:
            return [int(num_seq / n_chromosomes) + 10] * n_chromosomes
        except:
            print("Naive sampling needs the option --chromosomes...Aborting...", file=sys.stderr)
            sys.exit(1)


def main(infile, outfile, num_seq, length, chromosomes=None, min_length=35, max_length=100, unif=0.01):
    num_reads_to_sample = estimate_read_distribution(infile, num_seq, chromosomes)
    with gzip.open(infile, "rt") as file_in:
        all_chunks = []

        for num_record, record in enumerate(SeqIO.parse(file_in, "fasta")):
            # we want reads only to assigned chromosomes
            if chromosomes:
                print(f"Parsing chromosome {num_record + 1}", end="\r", file=sys.stderr)
                if num_record >= chromosomes:
                    break

            all_chunks += chunk_fast(
                record, num_reads_to_sample[num_record], unif=unif, len_distrib=length,
                minlength=min_length, maxlength=max_length)

        try:
            if len(all_chunks) > num_seq:
                all_chunks = random.sample(all_chunks, num_seq)
        except ValueError:
            print(f"Returning only {len(all_chunks)} sequences for {record.id}", file=sys.stderr)

    # create a new header which includes the read pos, read length
    record_it = (SeqRecord.SeqRecord(
        record.seq, id="{}|{}_{}".format(record.id, pos, len(record)),
        description=" ".join(record.description.split(' ')[1:])) for record, pos in sort_recs(all_chunks))

    with open(outfile, 'w') as file_out:
        SeqIO.write(record_it, file_out, "fasta")

    return True
