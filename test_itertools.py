#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import itertools

DNA_Seq = "TCAAATG"

def get_random(s, n):
    sample_N = []   
    num = len(s)
    N_base_comb = list(itertools.combinations(range(num), n))

    for primer_s in N_base_comb:
        L_nucl = list(s)
        for i in primer_s:           
            L_nucl[i] = "N"
        sample_N.append(L_nucl)
    return sample_N

test_s = get_random(DNA_Seq, 2)
print (test_s)
     
def generate_err_prone(pre_seq):
    NAN ={'W':"AT",'S':"GC",'M':"AC",'K':"GT",'R':"AG",'Y':"CT",'B':"CGT",'D':"AGT",'H':"ACT",'V':"ACG",'N':"ATCG"}
    L = [[i] for i in pre_seq]
    for i in L:
        for x in i:
            if x in NAN.keys():
                i.remove(x)
                for y in list(NAN[x]):
                    i.append(y)

    Err_prone_seqs = [[]]
    for i in L:
        t = []
        for y in i:
            for a in Err_prone_seqs:
                t.append(a+[y])
        Err_prone_seqs = t
    return Err_prone_seqs

def create_seq(test_list):
    for i in test_list:
        Degenerate_bases = generate_err_prone(i)     
        for j in Degenerate_bases:
            yield "".join(j)

for i in create_seq(test_s):
    print (i)
