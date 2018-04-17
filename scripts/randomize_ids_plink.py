#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 08:56:50 2018

@author: hhuang2
"""

## randomize IDs
from pandas_plink import read_plink
import pandas as pd
import numpy as np
#import argparse

## I/O test
#chrom = 22
#plink_fp = "../Data/" +str(chrom)+'bed'
#output = '../Data/'

## 
#parser = argparse.ArgumentParser()
#parser.add_argument("-i", "input", type=str,
#                    help="input file name with full directory.")
#parser.add_argument("-o", "output", type=str,
#                    help="output directory")
#args = parser.parse_args()
#
#plink_fp = args.input
#output = args.output
###

## Private cloud

plink_fp = '/mnt/cloudbiodata_nfs_2/users/hhuang/GWASH/cohort'
output = '/mnt/cloudbiodata_nfs_2/users/hhuang/GWASH/recode/'

##

(bim, fam, bed) = read_plink(plink_fp)

num_IDs = fam.shape[0]

new_ID_list = {'oldFID':['']*num_IDs,
               'oldIID':['']*num_IDs,
               'FID': ['']*num_IDs,
               'IID': ['']*num_IDs}
new_ID_list['oldFID'] = fam['fid'].tolist()
new_ID_list['FID'] = fam['fid'].tolist()
new_ID_list['oldIID'] = fam['iid'].tolist()

np.random.seed(1)
CaseIDs = np.random.permutation(range(1000000, 10000000))

new_ID_list['IID'] = ['RAND'+str(CaseIDs[x]) for x in range(num_IDs)]

new_ID_list_df = pd.DataFrame.from_dict(new_ID_list)
new_ID_list_df = new_ID_list_df[['oldFID', 'oldIID', 'FID', 'IID']]
new_ID_list_df.to_csv(output+'recode_list_cohort.tsv', sep = '\t', index = False, header = False)