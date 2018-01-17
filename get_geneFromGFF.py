#!/usr/bin/python


import sys, errno
import re
GFF_file = open(sys.argv[1],'r')
GFFC = GFF_file.readlines()
GFF_file.close()


for i in GFFC:
    cols = i.rstrip().split("\t")

    if i.startswith("#"):
        continue
    else:
        geneFunc = cols[1]
        regionFeature = cols[2]

        if geneFunc == "RefSeq" and regionFeature == "gene" :
            meta_info = cols[8]

            IDstrings = re.search(r'.*GeneID:(\d+?)[,;].*',meta_info)
            short_name = re.search('.*Name=([_a-zA-Z0-9\-\.\(\)\']+?);[\s\S]*',meta_info)
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




