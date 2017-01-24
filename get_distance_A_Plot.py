#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse



def get_D_Then_P(LD_f):
    pair_LD_e = [ 0, 0.1, 0.25, 0.4, 0.6, 1 ]
    d= {"0-0.025":{} , "0.025-0.05":{} , "0.05-0.075":{} , "0.075-0.12":{} , "0.12-0.2":{},
         "0.2-0.5":{}, "0.5-1.5":{}, "1.5-3":{}, "3-5": {},  "5-10":{}  }
    for wd in d:
        for pl in pair_LD_e[:-1]:
            d[wd][pl] = 0

    
    f = open(LD_f, 'r')
    for line in f:
        if not line.startswith(" CHR"):
            linSpl = line.rstrip().split()
            dist = int(linSpl[4]) - int(linSpl[1])
            r2 = float(linSpl[6])
            for n,j in enumerate(pair_LD_e[:-1]):
                if r2 <= pair_LD_e[n+1] and r2> j:
                    if dist >= 0 and dist < 25000:
                        d["0-0.025"][j] += 1
                    elif dist >= 25000 and dist < 50000:
                        d["0.025-0.05"][j] += 1
                    elif dist >= 50000 and dist < 75000:
                        d["0.05-0.075"][j] += 1
                    elif dist >= 75000 and dist < 120000:
                        d["0.075-0.12"][j] += 1
                    elif dist >= 120000 and dist < 200000:
                        d["0.12-0.2"][j] += 1
                    elif dist >= 200000 and dist < 500000:
                        d["0.2-0.5"][j] += 1
                    elif dist >= 500000 and dist < 1500000:
                        d["0.5-1.5"][j] += 1
                    elif dist >= 1500000 and dist < 3000000:
                        d["1.5-3"][j] += 1
                    elif dist >= 3000000 and dist < 5000000:
                        d["3-5"][j] += 1
                    elif dist >= 5000000 and dist < 10000000:
                        d["5-10"][j] += 1
    
    return d, pair_LD_e[:-1]

if __name__ == "__main__":

    USAGE = "Define input file and Total Sequencing Depth of each variant site"

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the refGene file containing SNP variations")
    parser.add_argument("-o", "--output", action = "store", help = "resultfile ")


    args = parser.parse_args()

    dr, PL = get_D_Then_P(args.input)
    
    if args.output != None:
        oh = open(args.output, 'w')
        oh.write("pair_wise_Distance")
        for i in sorted(PL, reverse=True):
            oh.write("\tA"+str(i))
        oh.write("\n")
    
    SortX = sorted(dr.items(), key= lambda x: float(x[0].split("-")[0])) 
    for i in SortX:
        if args.output != None:
            oh.write(i[0])
        else:
            print i[0]
        ks = sorted(i[1].keys(), reverse=True)
        for j in ks:
            #for k in range(dr[i[0]][j]):
                #print i[0], "A"+str(j)
            if args.output != None:
                oh.write(" {0}".format(dr[i[0]][j]))
            else:
                print dr[i[0]][j]
        if args.output != None:
            oh.write("\n")

    if args.output != None:
        oh.close()
    
    

