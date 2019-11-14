#!/usr/bin/env python
# -*- coding: UTF-8 -*-

## grep "V600E" *.xls | awk -F '\t' '{print $10"\t"$11"\t"$18}'
## awk -F '\t' '{print ($6"\t"$18"\t"$24)}' test1.txt

from collections import defaultdict
from glob import glob

f = defaultdict(dict)

a = open("mutations_celllines.txt", "r")
for line in a:
    x = line.rstrip().split("\t")
    if not line.startswith("Mutation_source"):
        f[x[0]][x[5]+"##" + x[3]] = float(x[8].replace("%", ""))

a.close()


f_test = defaultdict(dict)
mutFiles = glob("*_output.txt")
for mf in mutFiles:
    sample_name = mf.split("_")[0]
    ma = open(mf, "r")
    for line in ma:
        if not line.startswith("header"):
            mx = line.rstrip().split("\t")
            f_test[sample_name][mx[17] + "##" + mx[24]] = float(mx[5]) * 100
            #print(mx[17] + "##" + mx[23])

    ma.close()

for c in f:
    ref_dict = f[c]
    for p in ref_dict:

        if p in f_test[c]:
            genename = p.split("##")[0]
            AA_info  = p.split("##")[1]
            print("{0}\t{1}\t{2}\t{3}\tNGS_test {4}".format(c, genename, AA_info, ref_dict[p],f_test[c][p], ))


