#!/usr/bin/env python3 

gtf1709 = "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/info_files/hg19_genes_TE_1709.gtf"
gtf1811 = "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_syscol/info_files/hg19_genes_TE_1811_s.gtf"


dico = {}
with open(gtf1811,"rt") as f: 
    linecount = 0 
    for line in (line.rstrip().split("\t") for line in f):
        if linecount % 100000 == 0: print(f"Read {linecount}")
        linecount += 1
        geneID = line[8].split(";")[0].split(" ")[1].replace("\"","")
        key = f"{line[0]};{line[2]};{line[3]};{line[4]};{geneID}"
        dico[key] = ""

f.close()

Nsimilarities = 0 
Ndifferences = 0
g = open("differences1709_1811.txt","w")
with open(gtf1709,"rt") as f: 
    for line in (line.rstrip().split("\t") for line in f):
        geneID = line[8].split(";")[0].split(" ")[1].replace("\"","")
        key = f"{line[0]};{line[2]};{line[3]};{line[4]};{geneID}"
        #print(key)
        if key not in dico: 
            g.write("\t".join(line)+ "\n")
            Ndifferences += 1
        else: 
            Nsimilarities += 1
            
print(Ndifferences)
print(Nsimilarities)