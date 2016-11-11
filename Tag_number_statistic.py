#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
from collections import defaultdict
import argparse
import os

def UnFinished_samp(dir):
    files = os.listdir(dir)

    fin = []
    sample= []

    for i in files:
        if i.startswith("Unique_Sample_") and i.endswith("_res.fa"):
            fin.append(i.split("_")[2])
        elif i.startswith("sample_") and i.endswith(".fa"):
            sample.append(i.split("_")[1].rstrip(".fa"))
    print (sample)
    for i in sample:
        if i not in fin:
            retCod = os.popen("ls -l sample_%s.fa"%(i))
            print (retCod.read())

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
   
    USAGE = "python %s.py -i sum_file"%(__name__)
    parser = argparse.ArgumentParser(description=USAGE)     #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-i", "--input", action = 'store', dest='testName', required = True , help="Provide summary file containing tags number.")
    options = parser.parse_args()
    D = Tag_Num(options.testName)
    print ("Sample_Name\tTotal_Num\tAverage_Length")
    for k,i in D.items():
        print ("%s\t"%(k), end = "")
        num_L = ["", ""]
        for j in i:
            if j.startswith("Tot"):
                num_L[0] = str(round(float(j.split("\t")[1])))
            elif j.startswith("Ave"):
                num_L[1] = str(round(float(j.split("\t")[1])))
        print ("\t".join(num_L))
