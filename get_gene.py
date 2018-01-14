#!/usr/bin/python


import sys, errno
import re
GFF_file = open(sys.argv[1],'r')
GFFC = GFF_file.readlines()
GFF_file.close()

## purpose: generate bed file for genes gtf
## USAGE: python get_gene.py E_coli_K12.gtf > result.bed

## use grep command to check the gene number
## grep -n "\bexon\b" E_coli_K12.gtf | grep "protein_coding" | wc -l

for i in GFFC:
    cols = i.rstrip().split("\t")

    geneFunc = cols[1]
    regionFeature = cols[2]

    if not i.startswith("#") and geneFunc == "protein_coding" and regionFeature == "exon" :
        meta_info = cols[8]
#        col2 = i.rstrip().split("\t")[2]
#    print col3

        IDstrings = re.search(r'.*gene_id "(\w+\d+?)".*',meta_info)
        short_name = re.search(r'.*gene_name "([_a-zA-Z0-9\-\.\(\)]+?)"; .*',meta_info)

        geneID = IDstrings.group(1)
        geneName = short_name.group(1)

#       print (short_name.group(1))
        #bio_type = re.search(r'gene_biotype=(\w+)',meta_info)
        #if bio_type.group(1)=="protein_coding":
        #    try:
        #        Desc= re.search(r'.*description=(.+?);.*',meta_info).group(1)
        #    except:
        #        Desc="no_descri"
#        print (Desc)
       
        #bio = bio_type.group(1)
        strand = cols[6] 
        start = cols[3]
        end   = cols[4]
        
        try:
            print ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(cols[0],start, end, geneName, geneID, strand))


        except IOError as e:
            if e.errno == errno.EPIPE:
                sys.exit("line screening stopped.")
