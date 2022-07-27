#!/usr/bin/env python3

from sys import argv
import gzip
import tabix
import pandas as pd
import subprocess
from collections import defaultdict
import itertools 
import re 

#============================================================================================================================
DESC_COMMENT = "Script that extracts Snp dosages and gene/te expression and creates a file ready to be used for BN learning."
SCRIPT_NAME = "create_triplets.py"
#============================================================================================================================

"""
#===============================================================================
@author: Nikolaos Lykoskoufis
@date: 27 May 2020
@copyright: Copyright 2020, University of Geneva
Script that extracts Snp dosages and gene/te expression and creates a file ready to be used for BN learning.#===============================================================================
"""

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


class Readphenotypes(object):
    """
    [Class that reads phenotypes]
    """
    def __init__(self,bedFile):
        self.bedFile = bedFile

    def read_phenotypes(self):
        """[read phenotypes and creates dictionary containing for each feature CHR,TSS and strand information.]

        Returns:
            [nested dict] -- [dict of dict of feature information.]
        """        
        dico = defaultdict(dict)
        f = Utils.myopen(self.bedFile)
        header = f.readline()
        for line in (line.rstrip().split("\t") for line in f):
            dico[line[3]]["chr"] = line[0]
            dico[line[3]]["tss"] = line[2]
            dico[line[3]]["strand"] = line[5]
            dico[line[3]]['exp'] = line[6:]
        return dico

    def read_header(self):
        """
        [Return header of bedfile]

        Arguments:
            bedFile {string} -- [The bedfile to be read]

        Returns:
            [list] -- [header of the bed file]
        """        
        return Utils.myopen(self.bedFile).readline().rstrip().split("\t")


    def get_feature_inside_window(self,region):
        """[Extract all features inside region]

        Arguments:
            region {str} -- [chr:start-end]

        Returns:
            [list] -- [List of features found inside the region window.]
        """        
        command = f"tabix {self.bedFile} {region}"
        process = subprocess.check_output(command, shell=True).decode('utf-8')
        process = process.split("\n")[:-1]
        return process


class Vcf(object):

    def __init__(self,vcf_file):
        self.vcf_file = vcf_file 

    def read_header(self):
        """[read vcf and extract header using tabix from command line]

        Returns:
            [list] -- [a list containing the header of the VCF ##CHROM.]
        """        
        return subprocess.check_output("tabix -H {0} | grep '#CHROM'".format(self.vcf_file), shell=True).decode('utf-8').rstrip().split("\t")
    
    def extract_variant_inside_window(self, region,var_target):
        """[Extracts genotype information for specific variant]

        Arguments:
            region {[str]} -- ["1:1000-1001"]
            var_target {str} -- [rsid of snpID, whatever is given under the ID column in the VCF.]

        Returns:
            [list] -- [list containing genotype information about the variant in question.]
        """        
        url = self.vcf_file
        tb = tabix.open(url)
        records = tb.querys(region)
        records = [record for record in records]
        if len(records) >1 :
            return [i for i in records if i[2] == var_target][0]
        else :
            return records[0]


class CreateTriplets(object):

    def __init__(self,qtltools_result, vcf_file, bed_file):
        self.qtltools_result = qtltools_result
        self.vcf_file = vcf_file 
        self.bed_file = bed_file 
        
        #Initialize objects
        self.Bed = Readphenotypes(self.bed_file)
        self.VCF = Vcf(self.vcf_file)

        # Create dictionary for molecular phenotype expression
        print("\t* Reading phenotypes into dictionary")
        self.phen_dico = self.Bed.read_phenotypes()
        
        # Get samples
        print("\t* Getting Headers")
        self.bed_header = self.Bed.read_header()[6:]
        self.vcf_header = self.VCF.read_header()[9:]
      

    def create_triplets(self):      
        f = Utils.myopen(self.qtltools_result)
        for line in (line.rstrip().split(" ") for line in f):
            if line[0] == "phe_id":
                continue
            else:
                phe_id = line[0]
                gene = phe_id.split("_")[-1]
                te = "_".join(phe_id.split("_")[:-1])
    
                variant = line[7]
                variant_pos = f"{line[8]}:{line[9]}-{line[10]}"
                #print(variant_pos)
                # Get gene expression 
                gene_exp = self.phen_dico[gene]['exp']
                te_exp = self.phen_dico[te]['exp']

                # Get dosage                
                genotype = self.VCF.extract_variant_inside_window(variant_pos, variant)
                
                # Extracting dosage
                DS = [i.split(":")[1] for i in genotype[9:]]

                print(f"\t* Creating triplet file for {variant} {gene} {te}. For debugging purpose variant pos: {variant_pos}")

                outName = f"{variant}_{gene}_{te}_triplet_content.txt"
                
                df_exp = pd.DataFrame([gene_exp,te_exp,self.bed_header],index=[gene,te,"SampleID"])
                df_exp = df_exp.T
                df_DS = pd.DataFrame([self.vcf_header,DS],index=["SampleID",variant])
                df_DS = df_DS.T
                df = pd.merge(df_DS, df_exp, how="inner",on=['SampleID'])
                df.to_csv(outName,sep="\t",header=True,index=False)



#qtltools = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/OLIVIER_APPROACH/TUMOR/permute/permutation_te_gene_pairs_All_FDR5.significant.txt"

#vcf = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"

#/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EXPRESSION/TUMOR/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_filtered_50_percent_and_more_276_tumor_PC30_PC3_geno_normalTransform.bed.gz"

CreateTriplets(*argv[1:]).create_triplets()
