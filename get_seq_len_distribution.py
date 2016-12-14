#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from glob import glob
from collections import defaultdict
import numpy as np
import os
import argparse
import sys

def Finished_samp(fa_dir):
    if not os.path.isdir(fa_dir):
        sys.exit("Directory does not exist. Please check it.")
    files = os.listdir(fa_dir)
    fin = []
    sample= []

    for i in files:
        if i.startswith("Unique_Sample_") and i.endswith("_res.fa"):
            fin.append(i.split("_")[2])
        elif i.startswith("sample_") and i.endswith(".fa"):
            sample.append(i.split("_")[1].rstrip(".fa"))
    fa_list = ["Unique_Sample_" + i + "_res.fa" for i in fin]
    for i in sample:
        if i not in fin:
            retCod = os.popen("ls -l sample_%s.fa"%(i))
            print (retCod.read())

    return fa_list

def get_length_stat(f, stat_f):
    the_Seq_Rec = SeqIO.parse( f , "fasta")
    sample_name = ".".join(f.split(".")[:-1])
    len_Of_Seq = []
    Tot_Len = 0

    for record in the_Seq_Rec:
        length = len(record.seq)
        len_Of_Seq.append(length)
        Tot_Len += length
    arr = np.array(len_Of_Seq)

    mx_len = np.max(arr)
    mean_len = np.mean(arr)
    General_stat = open(stat_f, "a")
    General_stat.write("The Total number of Tags in sample {0} is: {1}\n".format(sample_name, arr.size))
    General_stat.write("The Min length of Tags in sample {0} is: {1}\n".format(sample_name, np.min(arr)))
    General_stat.write("The Max length of Tags in sample {0} is: {1}\n".format(sample_name, mx_len))
    General_stat.write("The Mean length of Tags in sample {0} is: {1}\n".format(sample_name, mean_len))
    General_stat.close()

    hist, bin_edges = np.histogram(arr, bins = np.arange(0, mx_len + 20, 20), density = False)
#print (hist)
#print (bin_edges)

    Matri_Len = {}

    for n,i in enumerate(bin_edges):
        try:
            Matri_Len[str(i) +"-" + str(i+20)] = hist[n]
            #print ("{0}-{1}\t{2}".format(i,i+20, hist[n]))
        except IndexError:
            pass
            #print ("{0}-{1}\texceed Number range".format(i, i+20))
    stat_result = open(sample_name + "_stat.txt", 'w')
    stat_result.write("The distribution of length is:\n")
    Sort_Mat = sorted(Matri_Len.iteritems(), key = lambda x: int(x[0].split("-")[0]))
    for n,v in Sort_Mat:
        stat_result.write("{0}\t{1}\n".format(n, v))
    stat_result.close()

def Tag_Num(summary_file):
    LEN_dis = defaultdict(list)
    try:
        File_Handle = open(summary_file, "r" )
        for line in File_Handle:
            lc = line.rstrip()
            if line.startswith("The Total number of Tags "):
                start_i = lc.index("Unique_Sample_")
                end_i = lc.index("_res ")
                S_Name_Tag = lc[start_i+14:end_i]
                LEN_dis[S_Name_Tag].append("Tot\t"+lc.split()[-1])
            elif line.startswith("The Mean length of Tags "):
                start_i = lc.index("Unique_Sample_")
                end_i = lc.index("_res ")
                S_Name_Tag = lc[start_i+14:end_i]
                LEN_dis[S_Name_Tag].append("Ave\t"+lc.split()[-1])
    except IOError:
        print ("Summary file not Found.")
    else:
        File_Handle.close()
    return LEN_dis



if __name__ == "__main__":
    parser = argparse.ArgumentParser()     #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-f", "--outputDir", type=str, action='store',  dest='FastaDir', default =".", help="input the dir of unique fasta files") 
    parser.add_argument("-a", "--Summ", type=str, action='store',  dest='summary_file', default ="General_all_s.txt", help="input of the summary file of all samples")
    parser.add_argument("-s", type=str, action='store',  dest='Summary_f',required = True, help="Provide a filename for the summary of tag numbers results")
    options = parser.parse_args()

    fas = Finished_samp(options.FastaDir) 

    if not os.path.exists(options.summary_file):
        for i in fas:
            sn = i[:-3]+"_stat.txt"
            if not os.path.exists(options.FastaDir + "/"+ sn):
                get_length_stat(options.FastaDir +"/"+i, options.summary_file)
    
    D = Tag_Num(options.summary_file)

    STs = open(options.Summary_f, 'w')   # STs is abbreviation of Sample_Tag_stat

    STs.write("Sample_Name\tTotal_Num\tAverage_Length\n")
    for k,i in D.items():
        STs.write("%s\t"%(k))
        num_L = ["", ""]
        for j in i:
            if j.startswith("Tot"):
                num_L[0] = str(round(float(j.split("\t")[1])))
            elif j.startswith("Ave"):
                num_L[1] = str(round(float(j.split("\t")[1])))
        STs.write("\t".join(num_L)+"\n")
    STs.close()



