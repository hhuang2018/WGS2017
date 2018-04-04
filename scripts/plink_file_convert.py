#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon April 2 22:16:54 2017

@author: hhuang2
"""

#import csv
from utils import coreFunctions as cf
from pandas_plink import read_plink
#from pandas import DataFrame
import numpy as np
#import graphviz
#from pandas_plink import example_file_prefix

plink_fp = "../Data/" 
metadata_fp = "../Data/available_cases.csv"
output = '../Data/'
#metadata_fp = "../Data/available_cases_HLA.csv"

chrList = range(1, 23) # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])
BMT_mm_table = np.array([], dtype = 'float64')
for chrom in chrList:
    try: 
        (bim, fam, bed) = read_plink(plink_fp+str(chrom)+'bed')
        # bim - pandas.DataFrame – Alleles.
        # fam - pandas.DataFrame – Samples.
        # G - numpy.ndarray – Genotype.

        # fam.iid.head()
        # G.head()
        for case_index in range(num_case):
            
            bmt_fams = fam.query("iid in ['" + metadata_avail_cases['NMDP_DID'][case_index] + 
                                           "', '" + 
                                           metadata_avail_cases['NMDP_RID'][case_index]+"']")
            gt = bed[bmt_fams.i.values, :].compute()
            bmt_gt = abs(gt[[0]] - gt[[1]])
            if BMT_mm_table.shape[0] == 0:
                BMT_mm_table = bmt_gt
            else:
                BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
        np.savez(output+'BMT_mm_table_chr_'+str(chrom)+'.npz', mm_table = BMT_mm_table, ID_table = metadata_avail_cases)
    except FileNotFoundError as e:
        print(e)



    
        
        
   