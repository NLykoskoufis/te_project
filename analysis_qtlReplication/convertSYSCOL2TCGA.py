#!/usr/bin/env python3 

import gzip 

dico = {}
with gzip.open("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_GTEx_TCGA_snpIDs.txt.gz", "rt") as f: 
    for line in (line.rstrip().split("\t") for line in f):
        dico[line[0]] = line[2]

f.close()

with open("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt", "rt") as f: 
    with open("TUMOR_teqtls.chrALL.significant.TCGAready.txt","w") as g: 
        g.write(f.readline())
        for line in (line.rstrip().split(" ") for line in f): 
            snpID = line[7]
            if dico.get(snpID) != None: 
                g.write(" ".join(line[:7]) + " "  + dico.get(snpID) + " " + " ".join(line[8:]) + "\n")


dico = {}
with gzip.open("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_GTEx_TCGA_snpIDs.txt.gz", "rt") as f:
    for line in (line.rstrip().split("\t") for line in f):
        dico[line[0]] = line[1]

f.close()



with open("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt", "rt") as f:
    with open("NORMAL_teqtls.chrALL.significant.GTExready.txt","w") as g:
        g.write(f.readline())
        for line in (line.rstrip().split(" ") for line in f):
            snpID = line[7]
            if dico.get(snpID) != None:
                g.write(" ".join(line[:7]) + " "  + dico.get(snpID) + " " + " ".join(line[8:]) + "\n")
