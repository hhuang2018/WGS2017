#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon April 2 22:16:54 2017

@author: hhuang2
"""

import sys

# sys.path.append('../utils/')
from utils import coreFunctions as cf
from pandas_plink import read_plink
import numpy as np
import pandas as pd
import zarr
import dask.array as da


## EC2
plink_fp = "/home/hhuang/efs/GWASH_IMPUTED_DATA/ImputeQC/"
metadata_fp = "/home/hhuang/data/GWAS/available_cases.csv"
output = '/home/hhuang/data/GWAS/BMT_PresAbs_mm/'

## EMR
#plink_fp = "/home/hadoop/data/"  # "s3://nmdp-plink-bucket/plink-files/ImputeQC/"
#metadata_fp = "/home/hadoop/data/available_cases.csv"  # "s3://nmdp-hhuang/GWAS_dev/metadata/available_cases.csv"
#output = '/home/hadoop/data/BMT_mm/'

# parser = argparse.ArgumentParser()

chrList = range(1, 23)  # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])

for chrom in chrList:
    BMT_mm_table_count = np.array([], dtype='float64')
    BMT_mm_table_count_dir = np.array([], dtype='float64')
    BMT_mm_table_presabs = np.array([], dtype='float64')

    try:
        (bim, fam, bed) = read_plink(plink_fp + str(chrom) + 'bed')
        # bim - pandas.DataFrame – Alleles.
        # fam - pandas.DataFrame – Samples.
        # G - numpy.ndarray – Genotype.

        # fam.iid.head()
        # G.head()
        for case_index in range(num_case):

            bmt_fams = fam.query("iid in ['" + metadata_avail_cases['NMDP_DID'][case_index] +
                                 "', '" +
                                 metadata_avail_cases['NMDP_RID'][case_index] + "']")
            gt = bed[:, bmt_fams.i.values].compute()
            # bmt_gt = abs(gt[[0]] - gt[[1]])
            bmt_gt_count = abs(gt[:, 0] - gt[:, 1]) # mismatch count based

            bmt_gt_count_dir = gt[:, 0] - gt[:, 1] # mismatch with direction

            bmt_gt_presAbs = bmt_gt_count
            bmt_gt_presAbs[bmt_gt_presAbs == 2] = 1 # mismatch presence-absence based

            if BMT_mm_table_count.shape[0] == 0:
                BMT_mm_table_count = bmt_gt_count
            else:
                # BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
                BMT_mm_table_count = np.vstack((BMT_mm_table_count, bmt_gt_count))

            BMT_mm_table_count_da = da.from_array(BMT_mm_table_count, chunks=1000)


            if BMT_mm_table_count_dir.shape[0] == 0:
                BMT_mm_table_count_dir = bmt_gt_count_dir
            else:
                # BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
                BMT_mm_table_count_dir = np.vstack((BMT_mm_table_count_dir, bmt_gt_count_dir))

            BMT_mm_table_count_dir_da = da.from_array(BMT_mm_table_count_dir, chunks=1000)

            if BMT_mm_table_presabs.shape[0] == 0:
                BMT_mm_table_presabs = bmt_gt_presAbs
            else:
                # BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
                BMT_mm_table_presabs = np.vstack((BMT_mm_table_presabs, bmt_gt_presAbs))

            BMT_mm_table_presabs_da = da.from_array(BMT_mm_table_presabs, chunks=1000)

        BMT_mm_table_count_da.to_zarr(output+'Encoded.zarr/'+str(chrom)+'/count')
        BMT_mm_table_count_dir_da.to_zarr(output+'Encoded.zarr/'+str(chrom)+'/count_directed')

        BMT_mm_table_presabs_da.to_zarr(output+'Encoded.zarr/'+str(chrom)+'/pres_abs')

    except FileNotFoundError as e:
        print(e)






