#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import math
import linecache
import argparse

def Dist_Matri(f, out):

    d = {}
    D_Ind = defaultdict(list)

    XT = 0
    geno_example = linecache.getline(f, 1)
    Split_geno = geno_example.rstrip().split(" ", 6)

    GeneType_List = Split_geno[6].split()
    GTL = []
    for n, t in enumerate(GeneType_List):
        if n%2 == 0:
            GTL.append(t)
            d[XT] = 0
            d["a" + str(XT)] = 0
            XT += 1
    
    ped_H = open(f, 'r')
    for i in ped_H:
        Sp = i.rstrip().split(" ", 6)
        ID = Sp[0]
        a = Sp[6]
        x = a.split()

        D_Ind[ID].append("")
        D_Ind[ID].append("")

        for num , base in enumerate(x):
            if num %2 == 0:
                loci = num / 2
                if base == GTL[loci]:
                    d[loci] += 1
                else:
                    d["a"+str(loci)] += 1
            else:
                loci = (num - 1) / 2
                if base == GTL[loci]:
                    d[loci] += 1
                else:
                    d["a"+str(loci)] += 1
    ped_H.close()

    N_ped = open(f, 'r')
    for i in N_ped:
        NSp = i.rstrip().split(" ", 6)
        ID = NSp[0]
        IN = NSp[6]
        xHap = IN.split()
        for num , base in enumerate(xHap):
            if num %2 == 0:
                Ind_loci = num / 2
                if base == GTL[Ind_loci]:
                    D_Ind[ID][0] += str(d[Ind_loci]) + "/"
                else:
                    D_Ind[ID][0] += str(d["a"+str(Ind_loci)]) + "/"
            else:
                Ind_loci = (num - 1) / 2
                if base == GTL[Ind_loci]:
                    D_Ind[ID][1] += str(d[Ind_loci]) + "/"
                else:
                    D_Ind[ID][1] += str(d["a"+str(Ind_loci)]) + "/"

    N_ped.close()

    N_G = D_Ind.keys()
    St_TSV_Handle = open(out, 'w')
    St_TSV_Handle.write("    {0}\n".format(len(N_G)))

    for num,ind in enumerate(N_G):
        #for n,j in enumerate(N_G[num+1:]):
        #St_TSV_Handle.write("{0}".format((ind, "<10")))
        St_TSV_Handle.write(str(ind).ljust(10))
        for n,j in enumerate(N_G):
            dist_L = 0
            X_g1 = D_Ind[ind][0].strip("/").split("/")
            X_g2 = D_Ind[ind][1].strip("/").split("/")
            Y_g1 = D_Ind[j][0].strip("/").split("/")
            Y_g2 = D_Ind[j][1].strip("/").split("/")
            for ni,i in enumerate(X_g1):
                X_a = int(X_g1[ni])
                X_b = int(X_g2[ni])
                Y_c = int(Y_g1[ni])
                Y_d = int(Y_g2[ni])

                Total = (d[ni] + d["a"+str(ni)]) * 1.0

                Delta_1 = math.sqrt((X_a/Total - Y_c/Total)**2) + math.sqrt((X_b/Total - Y_d/Total)**2)
                
                dist_L +=  Delta_1
            DistXY = "{0:.4f}".format(dist_L / len(GTL))
                
            St_TSV_Handle.write(DistXY.ljust(8))
        St_TSV_Handle.write("\n")
    St_TSV_Handle.close()
    return len(N_G)
        
    # ref  www.nature.com/articles/ncomms9824
    # ref  citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.478.3176&rep=rep1&type=pdf
    # ref  grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html

if __name__ == "__main__":

    USAGE = "python dist*.py -i inputFile -o outputFile"

    parser = argparse.ArgumentParser(description= USAGE)     #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-i", "--input", action = 'store', required=True, help="input a string passing to testName") # allows input of the file for reading
    parser.add_argument("-o", "--output", action = 'store', default="result_dist2.txt", help = "Output File name")
    options = parser.parse_args()

    a = Dist_Matri(options.input, options.output)


