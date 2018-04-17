#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 22:12:19 2018

@author: hhuang2
"""
import sys
sys.path.append('../utils/')
from utils import coreFunctions as cf
from pandas_plink import read_plink
#from pandas import DataFrame
import numpy as np

from pyspark import SparkContext as sc

#import argparse

#import graphviz
#from pandas_plink import example_file_prefix

#plink_fp = "../Data/" 
#metadata_fp = "../Data/available_cases.csv"
#output = '../Data/'

## EC2
#plink_fp = "/home/hhuang/efs/GWASH_IMPUTED_DATA/ImputeQC/"
#metadata_fp = "/home/hhuang/data/available_cases.csv"
#output = '/home/hhuang/data/BMT_mm/'

## EMR
#plink_fp = "/home/hadoop/data/" #"s3://nmdp-plink-bucket/plink-files/ImputeQC/"
#metadata_fp = "/home/hadoop/data/available_cases.csv"#"s3://nmdp-hhuang/GWAS_dev/metadata/available_cases.csv"
#output = '/home/hadoop/data/BMT_mm/'

#parser = argparse.ArgumentParser()

chrList = range(1, 23) # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])

for chrom in chrList:
    BMT_mm_table = np.array([], dtype = 'float64')
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
            gt = bed[:, bmt_fams.i.values].compute()
            #bmt_gt = abs(gt[[0]] - gt[[1]])
            bmt_gt = abs(gt[:, 0] - gt[:,1])
            if BMT_mm_table.shape[0] == 0:
                BMT_mm_table = bmt_gt
            else:
                # BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
                BMT_mm_table = np.vstack((BMT_mm_table, bmt_gt))
        np.savez(output+'BMT_mm_table_chr_'+str(chrom)+'.npz', mm_table = BMT_mm_table, ID_table = metadata_avail_cases)
    except FileNotFoundError as e:
        print(e)
        


    
        
        
   