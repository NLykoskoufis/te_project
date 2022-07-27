#!/usr/bin/env python3 

import gzip 
import re 

def readSNPconverter(ftxt,dataSet):
    dico = {}
    with gzip.open(ftxt,"rt") as f: 
        for line in (line.rstrip().split("\t") for line in f):
            key = line[1] if dataSet == "GTEx" else line[2]
            dico[key] = line[0]
    return dico 



def readBNresults(ftxt,fout,dico):
    with open(ftxt, "rt") as f: 
        with open(fout, "w") as g: 
            g.write(f.readline().rstrip() + "\tsyscol_triplet\n")
            for line in (line.rstrip().split("\t") for line in f):
                syscol_snpID = dico[line[0]]
                gene = line[1]
                te = re.sub("_chr.{1,2}_","_",line[2])
                syscol_triplet = f"{syscol_snpID};{gene};{te}"
                g.write("\t".join(line) + "\t" + syscol_triplet + "\n")
                
                
if __name__ == "__main__":
    dico = readSNPconverter("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_GTEx_TCGA_snpIDs.txt.gz","TCGA")
    readBNresults("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/tumor/SYSCOL_tumor_in_TCGA_COAD_teqtls.txt","/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/tumor/SYSCOL.TUMOR.in.TCGA_COAD.txt",dico)
    dico = readSNPconverter("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_GTEx_TCGA_snpIDs.txt.gz","GTEx")
    readBNresults("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/normal/SYSCOL_Normal_in_GTEx_ColonTransverse_teqtls.txt","/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/normal/SYSCOL.NORMAL.in.GTEx.ColonTransverse.txt",dico)



