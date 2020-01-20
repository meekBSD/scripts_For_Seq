

sampleFile = open( "23samples.txt",  "r")

out1 = open("specificity_results.xls", 'w')
out1.write("sample1\tsample2\tcfDNA_mutationNum\tffpe_mutNum\tall_MutNum\tsensitivity\tspec\taccu\tnegLik\tPPV\tNPV\n")

for line in sampleFile:
    if not line.startswith("血液编号"):

        k = line.rstrip().split("\t")

        if len(k) > 7:
            cfM= set(k[7].split())
            
            if len(k)>8:
            
                ffM= set(k[8].split())
            else:
                ffM = set([])

            allMut = list(cfM | ffM)
            
            var_test= []
            var_cf  = []

            for i in allMut:
                if i in ffM:
                    var_test.append(1)
                else:
                    var_test.append(0)
                if i in cfM:
                    var_cf.append(1)
                else:
                    var_cf.append(0)

            out1.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(k[0], k[1],len(cfM), len(ffM), len(allMut)))
            
            a= 0
            b = 0
            c  = 0
            d   = 0
            
            for i in allMut:
                if i in ffM and i in cfM:
                    a += 1
                if i in ffM and i not in cfM:
                    c += 1
                if i not in ffM and i in cfM:
                    b += 1
            
            if a +c >0:
                out1.write("敏感性\t{0}\t".format(a/(a+c)))
            else:
                out1.write("敏感性\tNA\t")
            if b >0:
                out1.write("特异性\t{0}\t".format(d/(b+d)))
            else:
                out1.write("特异性\tNA\t")
            out1.write("准确率\t{0}\t".format((a+d)/(a+b+c+d)))
            if b >0 and d > 0:
                
                out1.write("阴性似然比\t{0}\t".format((c/(a+c))/(d/(b+d))))
            else:
                out1.write("阴性似然比\tNA\t")
                
            out1.write("PPV\t{0}\t".format(a/(a+b)))
            
            if c > 0:
                
            
                out1.write("NPV\t{0}\n".format(d/(c+d)))
            
            else:
                out1.write("NPV\tNA\n")
            

        else:
            
            out1.write("{0}\t{1}\t-\t-\t-\n".format(k[0],k[1]))

sampleFile.close()
