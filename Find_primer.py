#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from Bio import SeqIO
from Bio.Seq import Seq
import os
import random
import argparse
import itertools

USAGE= "python %program -i fasta_file -e 3 -o test_out.fasta"

parser = argparse.ArgumentParser(description=USAGE) #simplifys the wording of using argparse as stated in the python tutorial
parser.add_argument("-i", type=str, action='store',  dest='input', help="input the fasta file for processing") # allows input of the fasta filename
parser.add_argument("-e", type=int, action='store',  dest='err_num', default = 2, help="input the false paired base number in primer")
#parser.add_argument("-l", type=int, action='store',  dest='Minimal_GeneLen', help="input the false paired base number in primer")
parser.add_argument("-F", type=str, action='store',  dest='Forw', default = "TTATTCATGCCGCCCTATCG", help="input the forward primer seq")
parser.add_argument("-R", type=str, action='store',  dest='Reve', default = "GCGGACATAGAACAGCTACT", help="input the reverse primer seq")
parser.add_argument("-o", "--output", help="Directs the output to a name of your choice")
args = parser.parse_args()


f = args.input
file_dest = open(args.output, 'w')

Forw = args.Forw
Reve = args.Reve

#Forw = "TTATTCATGCCGCCCTATCG"
#Reve = "GCGGACATAGAACAGCTACT"

log_file = open("process_log.txt",'w')

def get_reverse_com(s):
    basecom = {"A":"T","C":"G","G":"C","T":"A"}
    pre_letters = list(s)
    new_let = [basecom[base] for base in pre_letters]
    S_let = ''.join(new_let)
    return (S_let[::-1])

f_RC = get_reverse_com(Forw)
r_RC = get_reverse_com(Reve)


def get_random(s, n):
    sample_N = []
    
    num = len(s)
    N_base_comb = list(itertools.combinations(range(num), n))
    #primer_s = random.sample(list(range(num)), n)
    for primer_s in N_base_comb:
        L_nucl = list(s)
        for i in primer_s:
            
            L_nucl[i] = "N"
        sample_N.append(L_nucl)
    return sample_N
     
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


#print (test_L)

def create_seq(test_list):
    for i in test_list:
        Degenerate_bases = generate_err_prone(i)
        
        for j in Degenerate_bases:
            yield "".join(j)

for record in SeqIO.parse(f,'fasta'):
    seq = record.seq
    p_forward = get_random(Forw, args.err_num)[1]
    p_reverse = get_random(Reve, args.err_num)[1]
    not_primer = True 
    for i in create_seq(generate_err_prone(p_forward)):
        
        if seq.find(i) != -1:
            not_primer = False
            forw_ind = seq.find(i)
            if seq.find(r_RC) != -1:
                reve_ind = seq.find(r_RC)
                new_seq = seq[forw_ind: (reve_ind + len(Reve))]
                file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
            #elif forw_ind > args.Minimal_GeneLen + 100:
            elif forw_ind > 500:
            #    new_seq = seq[forw_ind: len(seq)]
            #    file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
                log_file.write("%s\t%s\n"%(">"+record.id,"Can't find reverse primer."))
            else:
            #    new_seq = seq[0: forw_ind + len(seq)]
            #    file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
                log_file.write("%s\t%s\n"%(">"+record.id,"Can't find reverse primer."))
            #print("%s\t%d\n%s"%(">"+record.id,len(new_seq),new_seq))
    for i in create_seq(generate_err_prone(p_reverse)):    
        if seq.find(i) != -1:
            not_primer = False
            r_RC_ind = seq.find(i)
            if seq.find(f_RC)!= -1:
                f_RC_ind = seq.find(f_RC)
                new_seq = seq[r_RC_ind: (f_RC_ind + len(Forw))]
                file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
            #elif f_RC_ind > args.Minimal_GeneLen + 100:
            elif f_RC_ind > 500:
            #    new_seq = seq[r_RC_ind: len(seq)]
            #    file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
                log_file.write("%s\t%s\n"%(">"+record.id,"Can't find forward primer."))
            else:
            #    new_seq = seq[0: r_RC_ind + len(Reve)]
            #    file_dest.write("%s\t%d\n%s\n"%(">"+record.id,len(new_seq),new_seq))
                log_file.write("%s\t%s\n"%(">"+record.id,"Can't find forward primer."))
            #print("%s\t%d\n%s"%(">"+record.id,len(new_seq),new_seq))

    if not not_primer == False:
        pass
        log_file.write("%s\t%s\n"%(">"+record.id,"Can't find any of primers."))

file_dest.close()
log_file.close()
