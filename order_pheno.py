#!/usr/bin/env python
# -*- coding:UTF-8 -*-

from glob import glob
from collections import defaultdict
import pandas as pd
import os
import argparse
import subprocess
import shlex
import shutil
import csv

USAGE = "python ord*.py -i inputFile -o outputFile"

parser = argparse.ArgumentParser(description= USAGE)     #simplifys the wording of using argparse as stated in the python tutorial

parser.add_argument("-i", "--input", action = 'store', default='my_pHap.ped', help="input a string passing to testName") # allows input of the file for reading
parser.add_argument("-p", "--pheno", action = 'store', required=True, help="input a string passing to testName") # allows input of the file for reading
parser.add_argument("-o", "--output", help = "Output File name", action= 'store')
parser.add_argument("-c", "--clust", action ='store', help = "designate a cov file for population stratification analysis")
options = parser.parse_args()

ped = options.input 
phe = options.pheno
out = options.output

def order_phen(ped, phe, out):
    ind_list = []
    x = open(ped, 'r')    # input file 1

    for i in x:
        iid = i.rstrip().split()[1]
        ind_list.append(iid)
    x.close()
    
    order_pheno = "Ordered_Phenotype.txt"    
    if out != None:
        result = open(out, 'w')   # output file
        pheno= open(phe, 'r').readlines()  # inputfile 2
        for i in pheno:
            if i.startswith("FID"):
                result.write(i)

        for i in ind_list:
            for j in pheno:
                if j.rstrip().split("\t")[1] == i:
                    result.write(j)
        result.close()
        order_pheno = out
    return ind_list, order_pheno

def P_F(f2):
    f_h = open(f2, 'r')
    geno_Dict = {}
    for Line in f_h:
        if Line[0:2]=="ID":
            Name=Line.strip()[3:]
            geno_Dict[Name] = ""
        else:
            try:
                #try to append genotype to individual ID with +=, assumes Name is a key
                geno_Dict[Name]+=Line
            except KeyError:    
                geno_Dict[Name]=Line
            except UnboundLocalError:
                continue
    f_h.close()
    return geno_Dict

def transpose(mat_file):
    with open(mat_file, 'r') as fin:
        rows = csv.reader(fin, delimiter = "\t", skipinitialspace = True)
        transposed = zip(*rows)
        with open("transposed.txt", "w") as fout:
            w = csv.writer(fout, delimiter="\t")
            w.writerows(transposed)

def trans_cov(cov_file):
    clus_D = {}
    cov_handle = open(cov_file, 'r')
    for i in cov_handle:
        cov_columns = i.rstrip().split()
        GN = cov_columns[2]     # GN is abbreviation of group name
        if cov_columns[2] not in clus_D:
            clus_D[GN] = cov_columns[1]
        else:
            clus_D[GN] += "\t" + cov_columns[1]
    cov_handle.close()
    return clus_D

def get_lines(file_obj):
    items = [""]*2
    for number, line in enumerate(file_obj):
        if line.startswith(" CHR"):
            SNP_line_n = number + 1
            items[0] = line.rstrip()

        elif items[0] != "" and number == SNP_line_n:
            items[1] = line.rstrip()
            yield (items)


def T_Structure(f, out):
    N_M = 0
    N_I = 0
    Marker_L = []
    Marker_H = open(f+".map", "r")
    for line in Marker_H:
        Marker_L.append(line.rstrip().split("\t")[1])
        N_M += 1
    Marker_H.close()
    geno_Int = {"A":"1", "C":"2", "G":"3", "T":"4", "?":"-1"}

    St_TSV_Handle = open(out + "structure_sample.txt", 'w')
    for i in Marker_L:
        St_TSV_Handle.write(" "+i)
    St_TSV_Handle.write("\n")

    ped_H = open(f+".ped", 'r')
    for i in ped_H:
        Sp = i.rstrip().split(" ", 6)
        N_I += 1
        a = Sp[6]
        if Sp[1] != Sp[0]:
            Pop = Sp[1]
        else:
            Pop = "1"
        x = a.split()

        A_geno = ""
        B_geno = ""
        for num , i in enumerate(x):
            if num %2 == 0:
                A_geno += geno_Int[i] + " "
            else:
                B_geno += geno_Int[i] + " "

        St_TSV_Handle.write("{0} {1} {2}\n".format(Sp[0], Pop, A_geno.rstrip()))
        St_TSV_Handle.write("{0} {1} {2}\n".format(Sp[0], Pop, B_geno.rstrip()))
    St_TSV_Handle.close()
    ped_H.close()
    return N_M, N_I

def MoveFile(source_F, dest_D):
    if not os.path.exists(dest_D):
        os.makedirs(dest_D)
    elif os.path.isfile(dest_D):
        os.makedirs(dest_D)
    else:
        pass

    Cur_Dir = os.getcwd()
    os.rename(source_F, os.path.join(Cur_Dir , dest_D, source_F))


def get_target_Lines(fr, s):

    f_handle = open(fr, "r")

    targets = []
    #flag = False
    for num, line in enumerate(f_handle):
        if line.startswith(s):
            tmp = num
        try:
            if num > tmp and line.rstrip() != "":
                targets.append(line.rstrip())
            elif num > tmp and line.rstrip() == "":
        #        flag = True
                break
        except NameError, e:
            pass

        #if flag:
        #    break

    return targets


# //////////////////////////////////////////////////////////////////////

ord_ids, phenofile = order_phen(ped, phe, out)
d = P_F("Myresult1_hapguess_switch.out")

genotype_num = 0

if not os.path.exists("new_phased.ped"):
    ZJ_ped = open("new_phased.ped", 'w')

    for i in ord_ids:
        diploid = d[i].strip().replace("?","0").split("\n")
        di_1 = diploid[0].split()
        di_2 = diploid[1].split()
        GTs = ""
        for num, base in enumerate(di_1):
            GTs += base + " " + di_2[num] + " "
        ZJ_ped.write("{0} {0} {1} {2} {3} {4} {5}\n".format(i, 0, 0 ,0, -9, GTs.strip(" ")))
        genotype_num += 1
    #print len(d[i])
    ZJ_ped.close()
    print ("Individuals' number is: {0}".format(genotype_num))

map_f = "new_phased.map"
if not os.path.exists(map_f):
    map_ori = ped.rstrip("ped")+ "map"
    map_file_handle = open(map_f, 'w')

    ori_handle = open(map_ori, 'r')
    SNP_ID = 0
    for i in ori_handle:
        ori_splits = i.rstrip().split("\t")
        if ori_splits[1] == ".":
            SNP_ID += 1
            map_file_handle.write("{0}\tXS{1}\t{2}\t{3}\n".format(ori_splits[0], SNP_ID, ori_splits[2], ori_splits[3]))
        else:
            map_file_handle.write(i)
    map_file_handle.close()
    #  awk -F '\t' 'BEGIN{ID=1}{if($2!="."){print $1"\t"$2"\t0\t"$4}else{print $1"\tNS"ID"\t0\t"$4;ID+=1}}' my_pHap.map > new_phased.map

    ori_handle.close()

# ../../plink --file new_phased --pheno ZA_new_order.txt --pheno-name Ratio2 --allow-no-sex --assoc --adjust --out GWAS_ratio2
plink = "/home/qzy/software/plink1.9/plink"

filter_f = "filtered_phased"

if not os.path.exists(filter_f + ".ped"):
    filter_cmd = "{0} --file new_phased --maf 0.02 --geno 0.05 --mind 0.05 --hwe 1e-3 --recode --out filtered_phased".format(plink)
    sf = shlex.split(filter_cmd)
    subprocess.call(sf)
if not os.path.exists("GWAS_filtered.bed"):
    GFB_CMD = shlex.split("{0} --file filtered_phased --make-bed --out GWAS_filtered".format(plink))
    subprocess.call(GFB_CMD)

phe_list = []
pheno_handle = open( phenofile , 'r')
for i in pheno_handle:
    if i.startswith("FID"):
        p = i.rstrip().split("\t")[2:]
        phe_list.extend(p)
pheno_handle.close()

cluster = options.clust   # previous use cov file

if cluster != None:
    if os.path.exists(cluster):
        trans_CD = trans_cov(cluster)
        if not os.path.exists(cluster+"_C.ped") or not os.path.exists(cluster+"_C.map"):
            newP_handle = open(cluster+"_C.ped", "w")
            example_ped = open(filter_f + ".ped", 'r')
            for ind in example_ped:
                ind_GL = ind.rstrip().split(" ", 6)           # ind_GL is abbreviation of individual Genotype line
                for i in trans_CD:
                    if ind_GL[1] in trans_CD[i].split("\t"):
                        first_Fcol = " ".join(ind_GL[:5])
                        newP_handle.write( "{0} {1} {2}\n".format(first_Fcol, i ,ind_GL[6]))
                        break
            newP_handle.close()
            example_ped.close()
            shutil.copy("filtered_phased.map", cluster+"_C.map")
        if not os.path.exists("gwas_Clr.bed"):
            cmd_bG = shlex.split("{0} --file {1} --make-bed --out gwas_Clr".format(plink, cluster+"_C" ))
            subprocess.call(cmd_bG)
        if not os.path.exists("ibd_matrix.mibs"):
            clr_m = shlex.split("{0} --bfile gwas_Clr --allow-no-sex --cluster --matrix --out ibd_{1}".format(plink, "matrix"))
            subprocess.call(clr_m)
        if not os.path.exists("clusterF.cluster2"):
        # --mh/--bd requires at least two valid clusters.
            clr_cmd = shlex.split("{0} --bfile gwas_Clr --allow-no-sex --cluster --mc 15 --ppc 0.05 --out clusterF".format(plink))
            subprocess.call(clr_cmd)
        if not os.path.exists("GWAS_cov_Clust.cmh"):
            asso_cmd = shlex.split("{0} --bfile gwas_Clr --allow-no-sex --mh --within {1} --adjust --out GWAS_cov_Clust".format(plink, "clusterF.cluster2"))
        # --mh and --mh2 require a case/control phenotype
            subprocess.call(asso_cmd)
    else:
        sys.exit("Cov file not found, please check its name.")
else:
    asso_cmd = "{0} --bfile test_bed --allow-no-sex --assoc --adjust --out ".format(plink)

    for i in phe_list:
        pheno_name = i.lower()
        to = open("transposed.txt", "r")
        for j in to:
            columns = j.rstrip().split("\t")
        #if columns[0] not in ["FID", "IID"]:
            if columns[0].lower() == pheno_name:
                orig_pedf = open(filter_f + ".ped", "r")
                phe_pedf = open(filter_f + "_" + pheno_name + ".ped", "w")
                for line_num, l in enumerate(orig_pedf):
                    ind_sp = l.rstrip().split(" ", 6)
                    sid = ind_sp[1]
                    phe_pedf.write("{0} {1} {2} {3} {4} {5} {6}\n".format(ind_sp[0], sid, 0, 0 ,0, columns[line_num+1], ind_sp[6]))
                phe_pedf.close()
                orig_pedf.close()
                shutil.copy("filtered_phased.map", filter_f + "_" + pheno_name + ".map")
                if not os.path.exists("gwas_" + pheno_name + ".bed"):
                    mk_bedCMD = shlex.split("{0} --file {1} --make-bed --out {2}".format(plink, filter_f + "_" + pheno_name, "gwas_"+pheno_name ))
                    m_pipe = subprocess.Popen(mk_bedCMD)
                    m_pipe.wait()
                c = shlex.split(asso_cmd.replace("test_bed", "gwas_" + pheno_name) + "GWAS_"+pheno_name)
    # "{0} --file {1} --pheno {2} --pheno-name {3}  --mh --within {4} --allow-no-sex --assoc --adjust --out GWAS_cov_{4}".format(plink, filter_f, phenofile,i, cov_file, pheno_name)
                pipe_plink = subprocess.Popen(c)
    # pipe_plink = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE)
                pipe_plink.wait()
            
        to.close()

if os.path.exists("ibd_matrix.mibs") and not os.path.exists("test_MDS.pdf"):
    R_MDS_Script = open("view_IBD_script.R", "w")
    R_MDS_Script.write("""#!/usr/bin/Rscript\n
       m <- as.matrix(read.table("ibd_matrix.mibs"))\n
       mds <- cmdscale(as.dist(1-m))\n
       k <- c(rep("green", 7), rep("blue", 168))\n
       pdf("test_MDS.pdf")\n
       plot(mds, pch=20, col=k)\n
       dev.off()\n

       """)
    R_MDS_Script.close()

    R_CMD = shlex.split("Rscript --slave --no-restore --no-save view_IBD_script.R")
    subprocess.call(R_CMD)

assoFiles = [i for i in os.listdir(".") if (i.endswith(".qassoc") or i.endswith(".cmh")) and not i.startswith("plot_data")]

for i in assoFiles:
    if not os.path.exists("Clump_results/") or (os.path.exists("Clump_results/") and os.listdir("Clump_results/") == []):
        clump_cmd = shlex.split("{0} --file filtered_phased --clump {1} --clump-verbose  --out {2}_clu".format(plink, i, i)) 
        subprocess.call(clump_cmd)
        with open(i + "_clu.clumped", 'r') as fobj, open("clump_sites_"+ i + ".txt", "w") as outClu:
            for site in get_lines(fobj):
            #print "\t".join(site)
                SNP_ID = site[1].split()[2]
                outClu.write(SNP_ID + "\n")

        Extract_SNP= shlex.split("{0} --bfile GWAS_filtered --extract clump_sites_{1}.txt --recode --out Selected_SNP_{1}".format(plink, i))
        subprocess.call(Extract_SNP)
        Marker_C, Ind_C = T_Structure("Selected_SNP_" + i, i)
        for kn in range(2,6):
            structure_CMD = shlex.split( "/home/qzy/software/console/structure -m /home/qzy/software/console/mainparams -e /home/qzy/software/console/extraparams -i {0}structure_sample.txt -K {1} -L {2} -N {3} -o {4}_outfile.txt".format(i, kn, Marker_C, Ind_C, i))
            subprocess.call(structure_CMD)
            tgets = get_target_Lines( i+"_outfile.txt_f"  , "Inferred ancestry")
            with open("plot_data" + i + "_" + str(kn), "w") as fo:
                fo.writelines([li+"\n" for li in tgets])
            
            if not os.path.exists("plot_Structure.R"):
                with open("plot_Structure.R", 'w') as psr:
                    psr.write('''#!/usr/bin/Rscript

            mt <- function(v){
                n <- length(v)
                t1 <- v[1:5]
                t2 <- as.numeric(v[6:n])
                tp <- sum(t2)

                if (tp > 1){
                    max_index <- which(t2==max(t2))[1]
                    t2[max_index] <- sprintf("%.3f", max(t2)-tp+1)
                } else if (tp < 1){
                    min_index <- which(t2==min(t2))[1]
                    t2[min_index] <- sprintf("%.3f", min(t2) + 1 - tp)
                }
                return(c(t1,t2))
            }

            pl_struct <- function(infile, out) {
                dat <- read.table(infile, skip = 1)
                x <- t(apply(dat, 1, mt)) 
                x1 <- as.data.frame(x) 
                data1 <- x1[order(x1[,6], decreasing=T),]
                pdf(paste(out,"pdf",sep="."), height=5, width=15)
                T <- ncol(data1)                                                          
                barplot(as.matrix(t(data1[,6:T])), names.arg=data1$V2, col=rainbow(4), border= rainbow(4), las=2, cex.names=0.35, cex.axis =0.65) 
                dev.off()
            }

            args <- commandArgs(TRUE)
            f1 <- args[1]
            f2 <- args[2]
            pl_struct(f1, f2)
                   ''')
            R_CMD_S = shlex.split("Rscript --slave --no-restore  plot_Structure.R plot_data{0}_{1} {0}_{1}_str".format(i, kn))
            subprocess.call(R_CMD_S)
            MoveFile(i+"_"+ str(kn) +"_str.pdf", "Structure_Results")
        
    if not os.path.exists(i + "_TOP_sig.txt"):
        d = pd.read_csv(i, delimiter= "\s+")
    # d = pd.read_csv("GWAS_water.qassoc", delimiter= "\s+")
        m_D = d.dropna()
        sort_D = m_D.sort_values(by="P").iloc[0:5]
        sort_C = sort_D.sort_values(by="BP")
    # sort_D.to_csv("result.txt", sep= " ")
        sort_C.to_csv(i +"_TOP_sig.txt", sep= " ", index=False)
    
SNP_files = glob("*TOP_sig.txt")

for i in SNP_files:
    SNP_dist = defaultdict(list)
    T_handle = open(i, "r")
    for line in T_handle:
        if not line.startswith("CHR"):
            allCols = line.rstrip().split()
            CHR = allCols[0] 
            Top_Sites = "\t".join(allCols[1:3])
            SNP_dist[CHR].append(Top_Sites)
    T_handle.close()
    print i + "#####################"

    for c in SNP_dist:
        SNP_IDs = [i.split("\t")[0] for i in SNP_dist[c]]
        BPs = [int(i.split("\t")[1]) for i in SNP_dist[c] ]
        tmp = BPs[0]
        print c+"\t"+SNP_IDs[0] + "\t" +  str(tmp)
        subprocess.call("{0} --bfile GWAS_filtered --snp {1} --window 250 --recode \
                        --out SNP_{1}".format(plink, SNP_IDs[0]), shell = True)
        for n, i in enumerate(BPs):
            if i - tmp >= 500000:
            
                print c+"\t"+SNP_IDs[n] + "\t" + str(i)
                subprocess.call("{0} --bfile GWAS_filtered --snp {1} --window 250 --recode \
                        --out SNP_{1}".format(plink, SNP_IDs[n]), shell = True)
                tmp = i
                
T_SNP_F = glob("SNP*.ped")
Haploview_bin = "/home/qzy/software/java_tools/Haploview.jar"

for i in T_SNP_F:
    File_Prefix = i.rstrip("*.ped")
    Hap_str = "java -cp /home/qzy/software/java_tools/ Haploview_ped_info -h {0} -p {1}.ped  -m {1}.map".format(Haploview_bin, File_Prefix)
    if not os.path.exists("LD_graph/") or (os.path.exists("LD_graph/") and os.listdir("LD_graph/") == []):
        Hap_cmd = shlex.split(Hap_str)
        subprocess.call(Hap_cmd)
        LD_stat = "new_" + i + ".LD"
        if os.path.exists(LD_stat):
            MoveFile(LD_stat, "LD_graph")
        if os.path.exists(LD_stat+ ".png"):
            MoveFile(LD_stat+ ".png", "LD_graph")

    if not os.path.exists("Block_results/") or (os.path.exists("Block_results/") and os.listdir("Block_results/") == []):
        SearchBlock = shlex.split("{0} --file {1} --blocks no-pheno-req --out block_{1}".format(plink , File_Prefix))
        subprocess.call(SearchBlock)
        if os.path.exists("block_"+ File_Prefix + ".blocks"):
            MoveFile("block_"+ File_Prefix + ".blocks", "Block_results")

    if not os.path.exists("LD_reports/") or (os.path.exists("LD_results/") and os.listdir("LD_results/") == []):
        LD_reportCMD = shlex.split("{0} --file {1} --r2 --ld-window-r2 0.5 --ld-window 5 --ld-window-kb 20 --out LD_Window20_{1}".format(plink, File_Prefix))
        subprocess.call(LD_reportCMD)
        if os.path.exists("LD_Window20_" + File_Prefix + ".ld"):
            MoveFile("LD_Window20_" + File_Prefix + ".ld", "LD_reports") 




