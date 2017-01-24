#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
from collections import defaultdict
import re
import time
import os

pattern_1 = '[\s\S]+mapped\s\((.*):-nan%\)'
pat_comp_1 = re.compile(pattern_1)

pattern_2 = '.*bam\sDepth:\s([\d|\.]+)'
pat_comp_2 = re.compile(pattern_2)

#test_rounds = 1000000

#begin_time = time.time()

#Dep_results = glob("*_Depth_Result.txt")
#Quali_stats = glob("qualimap_stat_*.txt")

def FindFile_sample(d_dir, pName):
    if not os.path.isdir(d_dir):
        sys.exit("Directory does not exist. Please check it.")
    T_files = []
    for root, dirs, files in os.walk(d_dir):
        for fr in files:
            # if fr.endswith("Depth_Result.txt"):
            if pName.find("*") != -1:
                name_s = pName.split("*")
                if fr.startswith(name_s[0]) and fr.endswith(name_s[1]):
                    f_path = os.path.join(root, fr)
                    T_files.append(f_path)
            elif fr.endswith(pName):
                f_path = os.path.join(root, fr)
                T_files.append(f_path)
    return T_files

C_D = os.getcwd()

All_files = FindFile_sample(C_D, "Depth_Result.txt")

Q_Fs = FindFile_sample(C_D, "qualimap_stat_*.txt")


d_sample = defaultdict(dict)

# for f in Dep_results:
for f in All_files:
    sample_name = f.split("/")[-1].split("_")[0]
    with open(f, 'r') as d_handle:
        for i in d_handle:
            a = pat_comp_1.match(i.rstrip())
            b = pat_comp_2.match(i.rstrip())
            if a != None:
                d_sample[sample_name]["Mapping Ratio"] = a.group(1)
                # print a.group(1)
            if b != None:
                d_sample[sample_name]["Average mapping Depth"] = b.group(1)
                # print b.group(1)

# for f in Quali_stats:
for f in Q_Fs:
    s_f = f.split("/")[-1].split("_")[-1].split(".")[0]
    with open(f, 'r') as Q_handle:
        for i in Q_handle:
            if i.startswith("Reference Size"):
                d_sample[s_f]["Reference Size"] = i.rstrip().split()[-1]
                # print i.rstrip()
            elif i.startswith("Region Size"):
                d_sample[s_f]["Region Size"] = i.rstrip().split()[-1]
                # print i.rstrip()
            elif i.startswith("Percentage of Reference"):
                d_sample[s_f]["Percentage of Reference"] = i.rstrip().split()[-1]
                # print i.rstrip()
            elif i.startswith("Mapped Reads Region"):
                MR_S = i.strip().lstrip("Mapped Reads Region").split("/")
                d_sample[s_f]["Mapped Region"] = MR_S[0].strip()
                d_sample[s_f]["Region mapping ratio"] = MR_S[1].strip()
                # print i.rstrip()
            elif i.startswith("Number of Reads"): 
                d_sample[s_f]["Number of Reads"] = i.rstrip().split()[-1]
                # print i.rstrip()

T_header = ["Reference Size","Number of Reads", "Mapping Ratio", "Region Size", "Percentage of Reference", "Mapped Region", "Region mapping ratio",
            "Average mapping Depth"]

result_DQ_handle = open("result_DQ.txt", 'w')

result_DQ_handle.write("{0}\t{1}\n".format("Sample_Name","\t".join(T_header)))
# sorted(d.items(), lambda x, y: cmp(x[1], y[1]), reverse=True) 
sortK = sorted(d_sample.keys())

regionS = d_sample[sortK[0]]["Region Size"]
Perc_Of_Ref = d_sample[sortK[0]]["Percentage of Reference"]

#for i in d_sample:
for i in sortK:
    #print i
    result_DQ_handle.write(i)
    for j in T_header:
        if j in d_sample[i]:
            # print d_sample[i][j]
            if j == "Region Size" and d_sample[i][j] != regionS:
            
                result_DQ_handle.write("\t{0}".format(regionS))
            elif j == "Percentage of Reference" and d_sample[i][j] != Perc_Of_Ref:
            
                result_DQ_handle.write("\t{0}".format(Perc_Of_Ref))
            else:
                result_DQ_handle.write("\t{0}".format(d_sample[i][j]))
    result_DQ_handle.write("\n")

result_DQ_handle.close()

#elapsed_time = time.time() - begin_time
#print(test_rounds * 2 / elapsed_time)
