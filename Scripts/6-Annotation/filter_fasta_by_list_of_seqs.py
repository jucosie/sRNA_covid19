#!/usr/bin/env python3
from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
seqs_set = set(line.strip() for line in open(sys.argv[2]))

for seq_record in ffile:
    try:
        seqs_set.remove(str(seq_record.seq))
    except KeyError:
        print(seq_record.format("fasta"), end='')
        continue
