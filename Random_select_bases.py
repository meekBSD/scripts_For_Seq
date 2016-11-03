#!/usr/bin/python

import random

def get_random(s, n):
    sample_N = []
    L_nucl = list(s)
    num = len(s)
    primer_s = random.sample(list(range(num)), n)
    for i in primer_s:
        sample_N.append(s[i])
        L_nucl[i] = "N"
    return (primer_s, L_nucl)

test_s = get_random("TTCCGTCGGGCATCGT", 2)[1]
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

test_L = generate_err_prone(test_s)
#print (test_L)

def create_seq(test_list):
    for i in test_list:
        yield "".join(i)
        
gen_seq = create_seq(test_L)
for i in gen_seq:
    print (i)
