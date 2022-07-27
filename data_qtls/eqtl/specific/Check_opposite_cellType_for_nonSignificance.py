#!/usr/bin/env python3

import gzip
from sys import argv
import pandas as pd
import glob
import re
import multiprocessing as mp
import argparse
import itertools
import numpy as np
from collections import defaultdict

# specific eQTL header results
#  gene	phe_chr	phe_start	phe_end	var_id	var_chr	var_start	var_end	cpm	cpm_nt	geno	dosage	tissue	ids	eQTL_pval	eqtl_int_beta	eqtl_normal_beta	eqtl_tumor_beta	eqtl_int_pval	eqtl_int_qval

#NORMAL_SPECIFIC_EQTL = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/NORMAL/SPECIFIC/normal_specific_eQTL_analysis_results_ALL_conditional.txt"
#TUMOR_SPECIFIC_EQTL  = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/TUMOR/SPECIFIC/tumor_specific_eQTL_analysis_results_ALL_conditional.txt"

#NORMAL_CONDITIONAL_EQTL = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/NORMAL/CONDITIONAL/conditional_chrALL_all_variants_normal.txt.gz"
#TUMOR_CONDITIONAL_EQTL = "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EQTL/CIS/TUMOR/CONDITIONAL/conditional_chrALL_all_variants.txt.gz"

#NORMAL_NOMINAL_EQTL = ""
#TUMOR_NOMINAL_EQTL = ""

class Utils:
    def __init__(self):
        pass

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
    
    def check_overlap(pos1,pos2,strd=False):
        """[Checking overlap of 0-based bed coordinates]

        Arguments:
            pos1 {[list]} -- [chr,start,end]
            pos2 {[list]} -- [chr,start,end]

        Keyword Arguments:
            strd {bool} -- [Whether to force strandeness or not] (default: {False})

        Returns:
            [str] -- [YES = overlap or NO = no overlap.]
        """
        chr1 = pos1[0]
        start1 = pos1[1]
        end1 = pos1[2]
            
        chr2 = pos2[0]
        start2 = pos2[1]
        end2 = pos2[2]
        
        strd1 = pos1[3]
        strd2 = pos2[3]
        # Checking whether we force strandeness or not. 
        if chr1 != chr2:
                raise Exception("You are trying to overlap genomic positions on different chromosomes. Are you sure there is no issue?")
        else:
            if strd == False:
                if start1 <= end2 and start2 <= end1: # if start of el1 is smaller than the end of el2 and the start of el2 is smaller than the end of el1 then there is an overlap. 
                    return "YES"
                else:
                    return "NO"
            else:
                if strd1 != strd2: # Even if they are overlapping, if the strand is different then no overlap. 
                    return "NO"
                else:
                    if start1 <= end2 and start2 <= end1:
                        return "YES"
                    else:
                        return "NO"
    
    def check_inside_phenotype(pos1,pos2,strd=False):
        """[Checking to see whether pos2 is inside pos1]

        Arguments:
            pos1 {[list]} -- [chr,start,end]
            pos2 {[list]} -- [chr,start,end]

        Keyword Arguments:
            strd {bool} -- [Whether to force strandeness or not] (default: {False})

        Returns:
            [str] -- [YES = inside or NO = outside.]
        """
        chr1 = pos1[0]
        start1 = pos1[1]
        end1 = pos1[2]
            
        chr2 = pos2[0]
        start2 = pos2[1]
        end2 = pos2[2]
        
        strd1 = pos1[3]
        strd2 = pos2[3]
        # Checking whether we force strandeness or not. 
        if chr1 != chr2:
            raise Exception("You are trying to overlap genomic positions on different chromosomes. Are you sure there is no issue?")
        else:
            if strd == False:
                if start2 >= start1 and end2 <= end1:
                    return "IN"
                else:
                    return "OUT"
            else:
                if strd1 != strd2:
                    return "OUT"
                else:
                    if start2 >= start1 and end2 <= end1:
                        return "IN"
                    else:
                        return "OUT"

class ReadQTL(object):

    def __init__(self,qtltools_output):
        self.qtltools_output = qtltools_output 
        
    def read_conditional_QTL(self):
        """
        [This function reads the conditional results from QTLtools cis and creates a dictionary for each eQTL_phen pair a dico with bwdpval.]
        """        
        dico = defaultdict(dict)
        f = Utils.myopen(self.qtltools_output)
        if f.readline().rstrip().split("\t")[0] == "phe_id":
            header = f.readline().rstrip().split("\t")
        for line in (line.rstrip().split(" ") for line in f):
            key = f"{line[7]};{line[0]}"
            # Populate dictionary 
            dico[key]['phe_id'] = line[0]
            dico[key]['phe_chr'] = line[1]
            dico[key]["phe_start"] = line[2]
            dico[key]['phe_end'] = line[3]
            dico[key]['phe_strand'] = line[4]
            dico[key]['n_var_in_cis'] = line[5]
            dico[key]['dist_phe_var'] = line[6]
            dico[key]['var_id'] = line[7]
            dico[key]['var_chr'] = line[8]
            dico[key]['var_start'] = line[9]
            dico[key]['var_end'] = line[10]
            dico[key]['rank'] = line[11]
            dico[key]['fwd_pval'] = float(line[12])
            dico[key]['fwd_r_squared'] = line[13]
            dico[key]['fwd_slope'] = float(line[14])
            dico[key]['fwd_best_hit'] = line[15]
            dico[key]['fwd_sig'] = line[16]
            dico[key]['bwd_pval'] = float(line[17])
            dico[key]['bwd_r_squared'] = float(line[18])
            dico[key]['bwd_slope'] = float(line[19])
            dico[key]['bwd_best_hit'] = line[20]
            dico[key]['bwd_sig'] = line[21]
        return dico

    def significant_eqtl_variants(self):
        dico = {}
        linecount = 0 
        f = Utils.myopen(self.qtltools_output)
        for line in (line.rstrip().split() for line in f):
            if len(line) > 1:
                raise Exception("There should be only one column containing the significant eQTL variants.")
            dico[line[0]] = ""
            linecount += 1
        print(f"Read {linecount} significant variants")
        return dico

    def nominal_eqtls(self):
        dico = defaultdict(dict)
        f = Utils.myopen(self.qtltools_output)
        for line in (line.rstrip().split(" ") for line in f):
            phe_id = line[0]
            var_id = line[7]
            nom_pval = line[11]
            slope = line[13]
            eqtl = f"{var_id};{phe_id}"
            dico[eqtl]["nom_pval"] = nom_pval 
            dico[eqtl][slope] = slope
        return dico 

    def nominal_eqtls_pd(self):
        df = pd.read_csv(self.qtltools_output, sep=" ",header=0)
        return df 
    
    def read_QTLspecific(self):
        f = Utils.myopen(self.qtltools_output)
        header = ""
        dico = defaultdict(dict)
        for line in (line.rstrip().split("\t") for line in f):
            if line[0] == "gene":
                header = line
            else: 
                phe_id = line[0]
                phe_chr = line[1]
                phe_start = line[2]
                phe_end = line[3]
                var_id = line[4]
                var_chr = line[5]
                var_start = line[6]
                var_end = line[7]
                cpm = line[8]
                cpm_nt = line[9]
                geno = line[10]
                dosage = line[11]
                tissue = line[12]
                ids = line[13]
                samples = line[14]
                eqtl_pval = line[15]
                eqtl_int_beta = line[16]
                eqtl_normal_beta = line[17]
                eqtl_tumor_beta = line[18]
                eqtl_int_pval = line[19]
                eqtl_int_qval = line[20]
                
                key = f"{var_id};{phe_id}"

                dico[key]['phe_id'] = phe_id 
                dico[key]['phe_start'] = phe_start 
                dico[key]['phe_chr'] = phe_chr 
                dico[key]['phe_end'] = phe_end 
                dico[key]['var_id'] = var_id 
                dico[key]['var_chr'] = var_chr
                dico[key]['var_start'] = var_start
                dico[key]['var_end'] = var_end 
                dico[key]['cpm'] = cpm 
                dico[key]['cpm_nt'] = cpm_nt
                dico[key]['geno'] = geno 
                dico[key]['dosage'] = dosage 
                dico[key]['tissue'] = tissue 
                dico[key]['ids'] = ids 
                dico[key]['samples'] = samples
                dico[key]['eqtl_pval'] = float(eqtl_pval)
                dico[key]['eqtl_int_beta'] = float(eqtl_int_beta)
                dico[key]['eqtl_normal_beta'] = float(eqtl_normal_beta)
                dico[key]['eqtl_tumor_beta'] = float(eqtl_tumor_beta)
                dico[key]['eqtl_int_pval'] = float(eqtl_int_pval)
                dico[key]['eqtl_int_qval'] = float(eqtl_int_qval)
        return dico 





# specific eQTL header results
#  gene	phe_chr	phe_start	phe_end	var_id	var_chr	var_start	var_end	cpm	cpm_nt	geno	dosage	tissue	ids	eQTL_pval	eqtl_int_beta	eqtl_normal_beta	eqtl_tumor_beta	eqtl_int_pval	eqtl_int_qval

def main(qtlSpecific_results, nominal_qtlFile, conditionalFile,output):
    """[File that filters for proper tissueSpecific eQTLs. If slopes are same sign, check to see whether the pvalue of association is not significant in the opposite tissue.]

    Args:
        qtlSpecific_results ([type]): [description]
        nominal_qtlFile ([type]): [description]
        conditionalFile ([type]): [description]
        output ([type]): [description]

    Returns:
        [None]: [The function write the significant results in a file.]
    """
    print(f"\t* Reading conditional results")
    tumor_conditional = ReadQTL(conditionalFile).read_conditional_QTL()
    print(f"\t* Reading tissue specific results")
    tissueSpe = ReadQTL(qtlSpecific_results).read_QTLspecific()

    g = open(output, "w")
    g.write("phe_id\tphe_chr\tphe_start\tphe_end\tphe_strand\tdummy\tdist_phe_var\tvar_id\tvar_chr\tvar_start\tvar_end\teqtl_int_beta\teqtl_normal_beta\teqtl_tumor_beta\tnominal_slope\tbwd_slope\tnominal_pval\tbwd_pval\teqtl_int_pval\teqtl_int_qval\tspe_hit\n")
    print(f"\t* Reading nominal results and filtering for proper tissue specific eQTLs")

    # Count variables 
    tissueSpe_significant = 0 
    tissueSpe_nonSignificant = 0 
    different_beta_sign = 0 
    same_beta_direction = 0
    same_beta_nominally_not_significant = 0 
    same_beta_nominally_significant = 0 

    f = Utils.myopen(nominal_qtlFile)
    line_count = 0 
    for line in (line.rstrip().split(" ") for line in f):
        if line[0] == "phe_id": continue
        line_count += 1 
        #print(f"{line_count}\r",end="")
        phe_id = line[0]
        var_id = line[7]
        nom_pval = float(line[11])
        nominal_slope = float(line[13])
        
        eqtl = f"{var_id};{phe_id}"

        # Check if nominal eQTL is in tissueSPe
        if eqtl not in tissueSpe:
            continue 
        else:
            if tissueSpe[eqtl]['eqtl_int_qval'] > 0.05:
                tissueSpe_nonSignificant += 1 # If pvalue is not significant continue
                ### IF not significant then spe_hit = 0
                
                g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n".format(tissueSpe[eqtl]['phe_id'],tissueSpe[eqtl]['phe_chr'],tissueSpe[eqtl]['phe_start'],tissueSpe[eqtl]['phe_end'],tumor_conditional[eqtl]['phe_strand'],".",tumor_conditional[eqtl]['dist_phe_var'],tissueSpe[eqtl]['var_id'],tissueSpe[eqtl]['var_chr'],tissueSpe[eqtl]['var_start'],tissueSpe[eqtl]['var_end'],tissueSpe[eqtl]['eqtl_int_beta'],tissueSpe[eqtl]['eqtl_normal_beta'],tissueSpe[eqtl]['eqtl_tumor_beta'],nominal_slope,tumor_conditional[eqtl]['bwd_slope'],nom_pval,tumor_conditional[eqtl]['bwd_pval'],tissueSpe[eqtl]['eqtl_int_pval'],tissueSpe[eqtl]['eqtl_int_qval'],"0"))

                
            else: # if interaction pvalue is significant, then:
                tissueSpe_significant += 1

                tumor_beta = tumor_conditional[eqtl]['bwd_slope']
                if np.sign(tumor_beta) != np.sign(nominal_slope): # if slopes are different sign
                    different_beta_sign += 1

                    ### IF significant and same sign then spe_hit = 1
                    g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n".format(tissueSpe[eqtl]['phe_id'],tissueSpe[eqtl]['phe_chr'],tissueSpe[eqtl]['phe_start'],tissueSpe[eqtl]['phe_end'],tumor_conditional[eqtl]['phe_strand'],".",tumor_conditional[eqtl]['dist_phe_var'],tissueSpe[eqtl]['var_id'],tissueSpe[eqtl]['var_chr'],tissueSpe[eqtl]['var_start'],tissueSpe[eqtl]['var_end'],tissueSpe[eqtl]['eqtl_int_beta'],tissueSpe[eqtl]['eqtl_normal_beta'],tissueSpe[eqtl]['eqtl_tumor_beta'],nominal_slope,tumor_conditional[eqtl]['bwd_slope'],nom_pval,tumor_conditional[eqtl]['bwd_pval'],tissueSpe[eqtl]['eqtl_int_pval'],tissueSpe[eqtl]['eqtl_int_qval'],"1"))
                else: # If slopes are same sign
                    if nom_pval < 0.05:
                        same_beta_nominally_significant += 1

                        ### IF significant can slopes same direction and nominal significant then spe_hit = 1
                        g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n".format(tissueSpe[eqtl]['phe_id'],tissueSpe[eqtl]['phe_chr'],tissueSpe[eqtl]['phe_start'],tissueSpe[eqtl]['phe_end'],tumor_conditional[eqtl]['phe_strand'],".",tumor_conditional[eqtl]['dist_phe_var'],tissueSpe[eqtl]['var_id'],tissueSpe[eqtl]['var_chr'],tissueSpe[eqtl]['var_start'],tissueSpe[eqtl]['var_end'],tissueSpe[eqtl]['eqtl_int_beta'],tissueSpe[eqtl]['eqtl_normal_beta'],tissueSpe[eqtl]['eqtl_tumor_beta'],nominal_slope,tumor_conditional[eqtl]['bwd_slope'],nom_pval,tumor_conditional[eqtl]['bwd_pval'],tissueSpe[eqtl]['eqtl_int_pval'],tissueSpe[eqtl]['eqtl_int_qval'],"0"))
                        
                    else:

                        same_beta_nominally_not_significant += 1 

                        ### IF significant and slopes same direction and nominal not significant then spe_hit = 1

                        g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n".format(tissueSpe[eqtl]['phe_id'],tissueSpe[eqtl]['phe_chr'],tissueSpe[eqtl]['phe_start'],tissueSpe[eqtl]['phe_end'],tumor_conditional[eqtl]['phe_strand'],".",tumor_conditional[eqtl]['dist_phe_var'],tissueSpe[eqtl]['var_id'],tissueSpe[eqtl]['var_chr'],tissueSpe[eqtl]['var_start'],tissueSpe[eqtl]['var_end'],tissueSpe[eqtl]['eqtl_int_beta'],tissueSpe[eqtl]['eqtl_normal_beta'],tissueSpe[eqtl]['eqtl_tumor_beta'],nominal_slope,tumor_conditional[eqtl]['bwd_slope'],nom_pval,tumor_conditional[eqtl]["bwd_pval"],tissueSpe[eqtl]['eqtl_int_pval'],tissueSpe[eqtl]['eqtl_int_qval'],"1"))

    print(f"\t* non significant tissue specific eQTLs: {tissueSpe_nonSignificant}")
    print(f"\t* Significant tissue specific eQTLs: {tissueSpe_significant}")
    print(f"\t\t* Opposite slope direction: {different_beta_sign}")
    print(f"\t\t* same slope direction significant in opposite condition: {same_beta_nominally_significant}")
    print(f"\t\t* same slope direction not significant in opposite condition: {same_beta_nominally_not_significant}")

    print(f"\t* Read {line_count} lines")
    return None 

main(*argv[1:])
