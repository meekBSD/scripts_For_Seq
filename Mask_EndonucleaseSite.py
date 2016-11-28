#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
import sys
import re

input_F= sys.argv[1]
output_F = sys.argv[2]

records = SeqIO.parse(input_F, 'fasta')
result = open(output_F, 'w')

'''
for i in records:
    result.write(">"+str(i.id)+"\n")
    seq_str = str(i.seq).replace("GTAC", "NNNN")
    result.write(seq_str+"\n")
'''

# finde GATC***** and replace this pattern as NNNNNNNNN

pattern = re.compile("(GTAC){1}[ATGC]{5}", flags=re.I)

for i in records:
    result.write(">" + str(i.id) + "\n")
    seq_str = re.sub(pattern, "NNNNNNNNN", str(i.seq))
    result.write(seq_str+"\n")

result.close()
