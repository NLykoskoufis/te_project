#!/usr/bin/env python3 

from sys import argv
import gzip
import tabix
import pandas as pd
import subprocess
from collections import defaultdict
import itertools 
import re 
import LDpairs_OOP
import argparse

#============================================================================================================================
DESC_COMMENT = "Script that finds all shared triplets between normal and tumor either with same eQTl variant or variants in LD. For the OLIVIER_APPROACH METHOD."
SCRIPT_NAME = "find_shared_triplets.py"
#============================================================================================================================

"""
#===============================================================================
@author: Nikolaos Lykoskoufis
@date: 14 May 2020
@copyright: Copyright 2020, University of Geneva
Script that finds all shared triplets between normal and tumor either with same eQTl variant or variants in LD.
#===============================================================================
"""

### STEPS ###

# 1. Read normal results in memory and create a dictionary containing as key the te_gene --> all information about triplet results. 
# 2. Loop over the tumor bnresults and find the common triplets. They are common IF:
#    2.1. The SNP in normal and tumor triplet is the same.
#    2.2. The normal SNP and tumor SNP are in high LD (Rsquared >= 0.90)
# 3. Write results in file and create dictionary containing tumor and normal results for the common triplet for later use.
# 4. Read the dictionary and get information about model shifts between normal and tumor as well as model percentages (Can be done in R as well.)


####### GLOBAL VARIABLES FOR DEBUGING #######

#NORMAL_EQTL = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/OLIVIER_APPROACH/NORMAL/permute/permutation_te_gene_pairs_All_FDR5.significant.txt"
#TUMOR_EQTL  = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/OLIVIER_APPROACH/TUMOR/permute/permutation_te_gene_pairs_All_FDR5.significant.txt"

#NORMAL_BN_RESULTS  = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/OLIVIER_APPROACH/NORMAL/bn/bnlearn_All_normal.txt.gz"
#TUMOR_BN_RESULTS   = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/OLIVIER_APPROACH/TUMOR/bn/bnlearn_All_tumor.txt.gz"

#VCF                = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"

#LD_THRESHOLD       = 0.90

#NORMAL_EQTL = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/NORMAL/CONDITIONAL/conditional_chrALL_top_variants_per_rank_normal.txt.gz"
#TUMOR_EQTL  = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/TUMOR/CONDITIONAL/conditional_chrALL_top_variants_per_rank.txt.gz"

#NORMAL_BN_RESULTS  = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/NORMAL/bn/sameLeadSNP_te_snp_bnlearn_All_processed.txt.gz"
#TUMOR_BN_RESULTS   = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/TUMOR/bn/sameLeadSNP_te_snp_bnlearn_All_tumor_processed.txt.gz"

#VCF                = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"

#LD_THRESHOLD       = 0.90


#############################################

class Utils:
    def __init__(self):
        pass
    
    @staticmethod
    def myopen(filename):
        """[checks to see whether file is compressed or not]

        Arguments:
            self {[string]} -- [file name]

        Returns:
            [object] -- [IO object]
        """        
        f = gzip.open(filename)
        try:
            f.read(2)
            f.close()
            return gzip.open(filename,"rt")
        except:
            f.close()
            return open(filename,"rt")
    
    @staticmethod
    def create_pairs(feature):
        """[Creates pairs between two features keeping always one gene in the pair.]

        Arguments:
            feature {list} -- [list of features]

        Returns:
            [list of lists] -- [list containing all pairs]
        """        
        r = re.compile('''ENSG''')
        comb = list(itertools.combinations(feature,2))
        pairs = []
        for i in comb:
            if any(r.match(f) for f in i) and len(re.findall(r,','.join(i)))<2:
                pairs.append(i)
        return pairs

    @staticmethod
    def getLD(snp1,pos1,snp2,pos2,vcf_file):
        snp1 = LDpairs_OOP.Snp(vcf_file,pos1,snp1)
        snp2 = LDpairs_OOP.Snp(vcf_file,pos2, snp2)
        LD = LDpairs_OOP.LD_calculation(vcf_file,snp1,snp2)
        return LD



class ReadBnlearnResults(object):

    def __init__(self,bnlearn_file):
        self.bnlearn_file = bnlearn_file
    

    @staticmethod 
    def getBest(model_probabilities):
        """[Gives the model with the highest probability between the three.]

        Args:
            model_probabilities ([list]): [List containing the probabilities of each of the three models [Causal, Reactive, Independent]]
        """        
        if model_probabilities.index(max(model_probabilities)) == 0:
            return "m1"
        elif model_probabilities.index(max(model_probabilities)) == 1:
            return "m2"
        else:
            return "m3"


    def readBnlearn(self):
        """[This function reads the results obtained from bnlearn.]

        Returns:
            [type]: [description]
        """        
        f = Utils.myopen(self.bnlearn_file)
        dico = defaultdict(dict)
        for line in (line.rstrip().split("\t") for line in f):
            if line[0] == "var_id" or line[0] == "variant": continue
            variant = line[0]
            gene = line[1]
            te = line[2]
            m1 = float(line[7])
            m2 = float(line[8])
            m3 = float(line[9])
            l1 = float(line[4])
            l2 = float(line[5])
            l3 = float(line[6])
            proba = [m1,m2,m3]
            
            #best = ReadBnlearnResults.getBest(proba)
            best = line[11]
            
            key = f"{gene}/{te}"
            dico[key]["variant"] = variant
            dico[key]["gene"] = gene 
            dico[key]["te"] = te
            dico[key]["proba"] = proba
            dico[key]["loglik"] = [l1,l2,l3]
            dico[key]["best"] = best 
            dico[key]["triplet"] = f"{variant}_{gene}_{te}"

        return dico
 
class ReadQTLtools(object):

    def __init__(self,qtltools_output):
        self.qtltools_output = qtltools_output 
    
                
    def read_permutation_variant(self):
            """
            [This function reads the conditional results from QTLtools cis and creates a dictionary with genomic position of the variant.]
            """        
            dico = {}
            f = Utils.myopen(self.qtltools_output)
            for line in (line.rstrip().split(" ") for line in f):
                if line[0] == "phe_id":
                    continue 
                var_id = line[7]
                var_pos = f"{line[8]}:{line[9]}-{line[10]}"
                dico[var_id] = var_pos 
            return dico 



def main(normal_bnresults, tumor_bnresults, normal_eqtl, tumor_eqtl, vcf_file, output_file, LD_THRESHOLD=0.9):
    print("\t* Reading normal bnlearn results")
    normal_bn = ReadBnlearnResults(normal_bnresults).readBnlearn()
    print("\t* Reading tumor qtls")
    teqtl = ReadQTLtools(tumor_eqtl).read_permutation_variant()
    print("\n* Reading normal qtls")
    neqtl = ReadQTLtools(normal_eqtl).read_permutation_variant()
    
    count_same_triplets = 0
    count_potential_common = 0
    count_LD_triplets = 0 


    tumor_causal = 0 
    tumor_reactive = 0 
    tumor_independent = 0 

    normal_causal = 0 
    normal_reactive = 0 
    normal_independent = 0 

    shifts = defaultdict(int)

    print("Initiating search for common triplets between normal and tumor")
    g = open(output_file, "w")
    g.write("normal_variant\ttumor_variant\tgene\tte\tnormal_triplet\ttumor_triplet\tnormal_m1\tnormal_m2\tnormal_m3\ttumor_m1\ttumor_m2\ttumor_m3\tnormal_L1\tnormal_L2\tnormal_L3\ttumor_L1\ttumor_L2\ttumor_L3\tnormal_best\ttumor_best\tRsquared_LD\n")
    f = Utils.myopen(tumor_bnresults)
    for line in (line.rstrip().split("\t") for line in f):
        if line[0] == "var_id" or line[0] == "variant":
            continue 
        tumor_variant = line[0]
        gene = line[1] 
        te = line[2]
        m1 = float(line[7])
        m2 = float(line[8])
        m3 = float(line[9])
        l1 = float(line[4])
        l2 = float(line[5])
        l3 = float(line[6])
        tumor_proba = [m1,m2,m3]
        #tumor_best = ReadBnlearnResults.getBest(tumor_proba)
        tumor_best  = line[11]
        tumor_triplet = f"{tumor_variant}_{gene}_{te}"   

        

        tumor_key = f"{gene}/{te}"
        if normal_bn.get(tumor_key) != None:
            if normal_bn.get(tumor_key)["triplet"] == tumor_triplet: 
                if normal_bn.get(tumor_key)["variant"] == tumor_variant: 
                    count_same_triplets += 1
                    ld_rsquared = 1.0
                    
                    # Get normal bn results 
                    normal_variant = normal_bn.get(tumor_key)["variant"]
                    normal_gene = normal_bn.get(tumor_key)["gene"]
                    normal_te = normal_bn.get(tumor_key)["te"]
                    normal_proba = normal_bn.get(tumor_key)["proba"]
                    normal_loglik = normal_bn.get(tumor_key)["loglik"]
                    normal_best = normal_bn.get(tumor_key)["best"]
                    normal_triplet = normal_bn.get(tumor_key)["triplet"]
                    if normal_gene != gene and normal_te != te:
                        raise Exception("There is an issue. Gene or TE do not match!!")
                    else:
                        if normal_best == "m1": normal_causal += 1 
                        if normal_best == "m2": normal_reactive += 1 
                        if normal_best == "m3": normal_independent += 1
                        
                        if tumor_best == "m1": tumor_causal += 1 
                        if tumor_best == "m2": tumor_reactive += 1 
                        if tumor_best == "m3": tumor_independent += 1
                        
                        model_shift = f"{normal_best}->{tumor_best}"
                        shifts[model_shift] += 1

                        g.write(f"{normal_variant}\t{tumor_variant}\t{gene}\t{te}\t{normal_triplet}\t{tumor_triplet}\t{str(normal_proba[0])}\t{str(normal_proba[1])}\t{str(normal_proba[2])}\t{str(tumor_proba[0])}\t{str(tumor_proba[1])}\t{str(tumor_proba[2])}\t{str(normal_loglik[0])}\t{str(normal_loglik[1])}\t{str(normal_loglik[2])}\t{str(l1)}\t{str(l2)}\t{str(l3)}\t{normal_best}\t{tumor_best}\t{ld_rsquared}\n")
            else:
                if normal_bn.get(tumor_key)["gene"] != gene or normal_bn.get(tumor_key)["te"] != te: 
                    raise Exception("Variant is not the same but the gene or TE as well. There is something wrong here.") 
                else:
                    count_potential_common += 1

                    normal_variant = normal_bn.get(tumor_key)["variant"]
                    normal_gene = normal_bn.get(tumor_key)["gene"]
                    normal_te = normal_bn.get(tumor_key)["te"]
                    normal_proba = normal_bn.get(tumor_key)["proba"]
                    normal_loglik = normal_bn.get(tumor_key)["loglik"]
                    normal_best = normal_bn.get(tumor_key)["best"]
                    normal_triplet = normal_bn.get(tumor_key)["triplet"]

                    normal_variant = normal_bn.get(tumor_key)["variant"]                
                    normal_variant_pos = neqtl[normal_variant]

                    tumor_variant_pos = teqtl[tumor_variant]

                    ld = Utils.getLD(normal_variant,normal_variant_pos,tumor_variant,tumor_variant_pos,vcf_file)
                    print(f"* LD between {normal_variant} and {tumor_variant} for {gene} and {te}: {ld.r2()}")

                    ld_rsquared = ld.r2()
                    if ld_rsquared >= LD_THRESHOLD:
                        count_LD_triplets += 1
                        g.write(f"{normal_variant}\t{tumor_variant}\t{gene}\t{te}\t{normal_triplet}\t{tumor_triplet}\t{str(normal_proba[0])}\t{str(normal_proba[1])}\t{str(normal_proba[2])}\t{str(tumor_proba[0])}\t{str(tumor_proba[1])}\t{str(tumor_proba[2])}\t{str(normal_loglik[0])}\t{str(normal_loglik[1])}\t{str(normal_loglik[2])}\t{str(l1)}\t{str(l2)}\t{str(l3)}\t{normal_best}\t{tumor_best}\t{ld_rsquared}\n")

                        if normal_best == "m1": normal_causal += 1 
                        if normal_best == "m2": normal_reactive += 1 
                        if normal_best == "m3": normal_independent += 1
                        if tumor_best == "m1": tumor_causal += 1 
                        if tumor_best == "m2": tumor_reactive += 1 
                        if tumor_best == "m3": tumor_independent += 1
                        model_shift = f"{normal_best}->{tumor_best}"
                        shifts[model_shift] += 1
        else:
            continue

    print(f"*\t Common triplets          : {count_same_triplets}")
    print(f"*\t triplet to test for LD   : {count_potential_common}")
    print(f"*\t triplets in LD (R2 >=0.9): {count_LD_triplets}")
    
    print(f"#Model Frequency\n###Normal###\n\tCausal: {normal_causal}\n\tReactive: {normal_reactive}\n\tIndependent: {normal_independent}\n###Tumor###\n\tCausal: {tumor_causal}\n\tReactive: {tumor_reactive}\n\tIndependent: {tumor_independent}")

    for key, value in shifts.items():
        print(key,"\t",value)



parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='Find shared triplets between normal and tumor')
parser.add_argument('-t', '--tumor', dest='TUMOR_BN_RESULTS',type=str,help="BN results for tumor triplets")
parser.add_argument('-n', '--normal', dest='NORMAL_BN_RESULTS',type=str,help="BN results for normal triplets")
parser.add_argument('-teqtl', '--tumor-eqtl', dest='TUMOR_EQTL',type=str,help="tumor eQTL results")
parser.add_argument('-neqtl', '--normal-eqtl', dest='NORMAL_EQTL',type=str,help="normal eQTL results")
parser.add_argument('-vcf', '--vcf', dest='VCF',type=str,help="VCF file for LD estimations")
parser.add_argument('-out', '--output', dest='OUTPUT_FILE', type=str, help='Output file')

parser.add_argument('-r2', '--ld-rsquared', dest='LD_THRESHOLD', type=float, help='LD threshold for filtering variants.')
args = parser.parse_args()


if __name__ == "__main__":
    if args.TUMOR_BN_RESULTS == None or args.NORMAL_BN_RESULTS == None or args.TUMOR_EQTL == None or args.NORMAL_EQTL == None or args.VCF == None or args.LD_THRESHOLD == None or args.OUTPUT_FILE == None:
        raise Exception("One argument is missing. Please use find_shared_triplets.py --help to see which arguments need to be specified.")
    else:
        main(args.NORMAL_BN_RESULTS, args.TUMOR_BN_RESULTS, args.NORMAL_EQTL, args.TUMOR_EQTL, args.VCF, args.OUTPUT_FILE,args.LD_THRESHOLD)
