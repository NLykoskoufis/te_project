#!/usr/bin/env python3

#import tabix
import math
import subprocess
import sys
import re

#===========================================================================================================
DESC_COMMENT = "Script that calculates LD between pairs of variants from 1000genomes. Modified from LDlink"
SCRIPT_NAME = "CODer.py"
#===========================================================================================================

"""
#===============================================================================
@author: Nikolaos Lykoskoufis
@date: 13 March 2020
@copyright: Copyright 2020, University of Geneva
Script to calculate LD between pairs of variants from 1000genomes.
#===============================================================================
"""

class Snp(object):

    #### GLOBAL VARIABLES #####

    RSID_SYSCOL_FILE = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/rsid_snpID.bed.gz"


    def __init__(self, vcf_file, var_region, rsid):
        self.vcf_file = vcf_file
        self.var_region = var_region
        self.rsid = rsid
        self.GT = self.get_variant()


    def get_rsid(self):
        return self.rsid

    def get_pos(self):
        return self.var_region

    def syscol_rsid_to_snpID_1000g(FILE,region,snpID_1000g):
        '''
        This function reads syscol_rsid_file and checks to see whether the rsid in 1000genomes and syscol are pointing to the same variant.
        '''
        command = f"tabix {FILE} {region}"
        process = subprocess.check_output(command, shell=True).decode('utf-8')
        if len(process.split("\n")[:-1]) > 1:
            for line in (line.rstrip().split("\t") for line in process):
                if line[4] == snpID_1000g:
                    return True
        elif len(process.split("\n")) == 0:
            return False
        else:
            process = process.rstrip().split("\t")
            if process[4] == snpID_1000g:
                return True
            else:
                return False

    def get_variant(self):
        '''
        This function reads an indexed vcf_file or bcf_file and extract the variant inside the given region. It also checks that the rsid of the extracted variants is the same or if not that the snpID are similar.
        '''
        command = f"tabix {self.vcf_file} {self.var_region}"
        process = subprocess.check_output(command, shell=True).decode('utf-8')
        process = process.split("\n")[:-1]
        if len(process) > 1:
            for line in (line.rstrip().split("\t") for line in process):
                if line[2] == self.rsid: #or Snp.syscol_rsid_to_snpID_1000g(Snp.RSID_SYSCOL_FILE,self.var_region,f"{line[0]}_{line[1]}_{line[3]}_{line[4]}_hg37"):
                    return [i.split(":")[0] for i in line]
                else:
                    continue
        elif len(process) == 0:
            raise Exception("no variants1")
        else:
            # print(process)
            process = process[0].split("\t")
            if process[2] == self.rsid: # or Snp.syscol_rsid_to_snpID_1000g(Snp.RSID_SYSCOL_FILE,self.var_region,f"{process[0]}_{process[1]}_{process[3]}_{process[4]}_hg37"):
                return [i.split(":")[0] for i in process]
            else:
                raise Exception("Problem")

    def check_variant(var):
        '''
        Returns true if variants are biallelic, else false.
        '''
        if "," not in var[3] or "," not in var[4]:
            return True

    def prepare_alleles(self):
        '''
        This function extracts the genotypes and creates a dictionary that contains for all three genotypes the corresponding alleles.
        '''
        geno = self.GT
        if Snp.check_variant(geno) != True:
            raise Exception("Variant is not biallelic")
            sys.exit(1)
        if len(geno[3]) == 1 and len(geno[4]) == 1:
            snp1_a1 = geno[3]
            snp1_a2 = geno[4]
        elif len(geno[3]) == 1 and len(geno[4]) > 1:
            snp1_a1 = "-"
            snp1_a2 = geno[4][1:]
        elif len(geno[3]) > 1 and len(geno[4]) == 1:
            snp1_a1 = geno[3][1:]
            snp1_a2 = "-"
        elif len(geno[3]) > 1 and len(geno[4]) > 1:
            snp1_a1 = geno[3][1:]
            snp1_a2 = geno[4][1:]

        allele = {"0|0": [snp1_a1, snp1_a1], "0|1": [snp1_a1, snp1_a2], "1|0": [snp1_a2, snp1_a1], "1|1": [
            snp1_a2, snp1_a2], "0": [snp1_a1, "."], "1": [snp1_a2, "."], "./.": [".", "."], ".": [".", "."]}
        return allele

    def get_alleles(self):
        '''
        Extracts alleles for variant
        '''
        geno = self.GT
        if Snp.check_variant(geno) != True:
            raise Exception("Variant is not biallelic")
            sys.exit(1)
        if len(geno[3]) == 1 and len(geno[4]) == 1:
            snp1_a1 = geno[3]
            snp1_a2 = geno[4]
        elif len(geno[3]) == 1 and len(geno[4]) > 1:
            snp1_a1 = "-"
            snp1_a2 = geno[4][1:]
        elif len(geno[3]) > 1 and len(geno[4]) == 1:
            snp1_a1 = geno[3][1:]
            snp1_a2 = "-"
        elif len(geno[3]) > 1 and len(geno[4]) > 1:
            snp1_a1 = geno[3][1:]
            snp1_a2 = geno[4][1:]
        return [snp1_a1, snp1_a2]


class LD_calculation(object):

    #### GLOBAL VARIABLE ######
    POP_FILE = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/CAUSAL_INF/SYSCOL.list"
    re_sample = re.compile('''[NTP].{0,1}''')


    def __init__(self,vcf_file,snp1,snp2):
        self.vcf_file = vcf_file
        self.snp1 = snp1
        self.snp2 = snp2
        self.vcf_header = [LD_calculation.re_sample.sub("",i) for i in subprocess.check_output("tabix -H {0} | grep '#CHROM'".format(vcf_file), shell=True).decode('utf-8').rstrip().split("\t")]
        self.pop = LD_calculation.get_pop(LD_calculation.POP_FILE)

    def get_pop(pop_file):
        '''
        Later you will be able to choose from populations inside the 1000genomes. But for now and lack of time i am sticking with samples from SYSCOL.
        '''
        f = open(pop_file,"rt")
        pop = {}
        for line in (line.rstrip().split() for line in f):
            pop[line[0]] = ""
        return pop

    def combine_haplotypes(self):
        head = self.vcf_header
        alleles1 = self.snp1.prepare_alleles()
        genotype1 = self.snp1.get_variant()
        alleles2 = self.snp2.prepare_alleles()
        genotype2= self.snp2.get_variant()
        geno = {}
        for i in range(9, len(head)):
            geno[head[i]] = [alleles1[genotype1[i]],".."]

        for i in range(9,len(head)):
            if head[i] in geno :
                geno[head[i]][1] = alleles2[genotype2[i]]
        return geno

    def extract_haplotypes(self):

        pop_ids = self.pop
        geno = LD_calculation.combine_haplotypes(self)
        snp1_alleles = self.snp1.get_alleles()
        snp2_alleles = self.snp2.get_alleles()

        hap = {}
        statistics = {}

        for ind in pop_ids:
            if ind in geno:
                hap1 = geno[ind][0][0] + "_" + geno[ind][1][0]
                hap2 = geno[ind][0][1] + "_" + geno[ind][1][1]
                if hap1 in hap:
                    hap[hap1] += 1
                else:
                    hap[hap1] = 1

                if hap2 in hap:
                    hap[hap2] += 1
                else:
                    hap[hap2] = 1
        # Remove missing haplotypes
        keys = list(hap.keys())
        for key in keys:
            if "." in key:
                hap.pop(key, None)

        # Check all haplotypes are present
        if len(hap) != 4:
            snp1_a = [snp1_alleles[0], snp1_alleles[1]]  # [snp1_a1, snp1_a2]
            snp2_a = [snp2_alleles[0], snp2_alleles[1]]  # [snp2_a1, snp2_a2]
            haps = [snp1_a[0] + "_" + snp2_a[0], snp1_a[0] + "_" + snp2_a[1],
                    snp1_a[1] + "_" + snp2_a[0], snp1_a[1] + "_" + snp2_a[1]]
            for i in haps:
                if i not in hap:
                    hap[i] = 0

        # Sort haplotypes
        A = hap[sorted(hap)[0]]
        B = hap[sorted(hap)[1]]
        C = hap[sorted(hap)[2]]
        D = hap[sorted(hap)[3]]
        N = A + B + C + D
        tmax = max(A, B, C, D)

        hap1 = sorted(hap, key=hap.get, reverse=True)[0]
        hap2 = sorted(hap, key=hap.get, reverse=True)[1]
        hap3 = sorted(hap, key=hap.get, reverse=True)[2]
        hap4 = sorted(hap, key=hap.get, reverse=True)[3]

        delta = float(A * D - B * C)

        Ms = float((A + C) * (B + D) * (A + B) * (C + D))

        if Ms != 0:

            # D prime
            if delta < 0:
                D_prime = abs(delta / min((A + C) * (A + B), (B + D) * (C + D)))
            else:
                D_prime = abs(delta / min((A + C) * (C + D), (A + B) * (B + D)))

            # R2
            r2 = (delta ** 2) / Ms

            # P-value
            num = (A + B + C + D) * (A * D - B * C) ** 2
            denom = Ms
            chisq = num / denom
            p = 2 * (1 - (0.5 * (1 + math.erf(chisq ** 0.5 / 2 ** 0.5))))

        else:
            D_prime = "NA"
            r2 = "NA"
            chisq = "NA"
            p = "NA"

        statistics["D_prime"] = D_prime
        statistics["r2"] = r2
        statistics["chisq"] = chisq
        statistics["p"] = p

        return statistics

    def r2(self):
        res = LD_calculation.extract_haplotypes(self)
        return res["r2"]

    def D_prime(self):
        res = LD_calculation.extract_haplotypes(self)
        return res["D_prime"]

    def chisq(self):
        res = LD_calculation.extract_haplotypes(self)
        return res["chisq"]

    def pvalue(self):
        res = LD_calculation.extract_haplotypes(self)
        return res["p"]



def run(vcf_file,rsid1,region1,rsid2,region2):
    snp1 = Snp(vcf_file, region1, rsid1)
    snp2 = Snp(vcf_file, region2, rsid2)
    LD= LD_calculation(vcf_file,snp1,snp2)
    print(f"f\t{rsid1}\t{region1}\t{rsid2}\t{region2}\t{str(LD.r2())}\t{str(LD.pvalue())}\t{str(LD.chisq())}\t{str(LD.D_prime())}\t0")

if __name__ == "__main__":
    run(*sys.argv[1:])




#vcf_file = "/srv/beegfs/scratch/users/l/lykoskou/DATA/1000g/BCF/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf"
#vcf_file = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
#vcf_file = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
#region = "10:100163120-100163120"
#rsid = "rs10883089"

#region2= "10:1282828-1282828"
#rsid2= "rs3793735"

#snp1 = Snp(vcf_file, region, rsid)
#print(snp1.get_variant())
#snp2 = Snp(vcf_file, region2, rsid2)
#print(snp2.get_variant())
#LD= LD_calculation(vcf_file,snp1,snp2)
#print(LD.extract_haplotypes())
#print(LD.r2())
#rs2244500  18:661005-661005      rs9956056 18:5237162-5237162
#0.003943307128864339

#rs2244500  18:661005-661005      rs2846970 18:14133493-14133493
#0.0016070145439186071