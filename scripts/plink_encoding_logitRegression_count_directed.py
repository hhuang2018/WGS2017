#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur July 19 00:00:00 2018

@author: hhuang2
"""

import sys

# sys.path.append('/home/hhuang/WGS2018/utils/')
from utils import coreFunctions as cf
import numpy as np
import pandas as pd

from pandas_plink import read_plink

import gc
import rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

import pickle


# R packages
r = ro.r
stats = importr('stats')
base = importr('base')
rGlimnnet = importr('glmnet')
Vector = rpy2.robjects.vectors.IntVector
StrVector = rpy2.robjects.vectors.StrVector

# define R object function
r('''
        # create a function `LogisticRegression`
        LogitRegression <- function(x, y) {
            x[x==-1] = NA
            snpData = cbind(x, y)

            snpData = as.data.frame(snpData)
            colnames(snpData) = c('gt', 'Outcome')

            Logit_regression <- glm(formula = Outcome ~ gt, family = binomial(), data = snpData, na.action = na.exclude)

            p_val <- coef(summary(Logit_regression))[2,'Pr(>|z|)']
            return(p_val)
            # return(summary(regression))
        }
        # call the function `f` with argument value 3

        ''')

r_logitRegression = ro.globalenv['LogitRegression']

mode = 'count_directed' # count; pres_abs; count_directed

## EC2
plink_fp = "/home/hhuang/efs/GWASH_IMPUTED_DATA/ImputeQC/"
#Zarr_fp = "/home/hhuang/data/GWAS/BMT_PresAbs_mm/Encoded.zarr/"
metadata_fp = "/home/hhuang/data/GWAS/available_cases.csv"
output_fp = '/home/hhuang/data/GWAS/MismatchEncoded/' + mode + '/'

chrList = range(5, 23)  # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])

metadata_pd = pd.read_csv(metadata_fp, index_col=2)
metadata_pd.index = metadata_pd.index.map(str)
#agvhd24_labels = metadata_pd.loc[:,'agvhi24']
#agvhd34_labels = metadata_pd.loc[:,'agvhi34']
#cgvhd_labels = metadata_pd.loc[:,'cgvhi']


for chrom in chrList:
    # BMT_mm_table_count_dir = np.array([], dtype='float64')

    try:
        (bim, fam, bed) = read_plink(plink_fp + str(chrom) + 'bed')
        # bim - pandas.DataFrame – Alleles.
        # fam - pandas.DataFrame – Samples.
        # G - numpy.ndarray – Genotype.

        BMT_mm_table_count_dir = np.empty([num_case, bim.shape[0]], dtype='float64') - 1

        # fam.iid.head()
        # G.head()

        for case_index in range(num_case):
            bmt_fams_d = fam.query("iid in ['" + metadata_avail_cases.iloc[case_index]['did'] + "']")

            bmt_fams_r = fam.query("iid in ['" + metadata_avail_cases.iloc[case_index]['rid'] + "']")

            gt_d = bed[:, bmt_fams_d.i.values].compute()
            gt_r = bed[:, bmt_fams_r.i.values].compute()

            # gt = bed[:, bmt_fams.i.values].compute()
            # bmt_gt = abs(gt[[0]] - gt[[1]])

            # bmt_gt_count_dir = gt[:, 0] - gt[:, 1]  # mismatch with direction

            bmt_gt_count_dir = gt_d - gt_r  # mismatch with direction

            BMT_mm_table_count_dir[case_index, :] = bmt_gt_count_dir
            #if BMT_mm_table_count_dir.shape[0] == 0:
            #    BMT_mm_table_count_dir = bmt_gt_count_dir
            #else:
                # BMT_mm_table = np.concatenate((BMT_mm_table, bmt_gt), axis =0)
            #    BMT_mm_table_count_dir = np.vstack((BMT_mm_table_count_dir, bmt_gt_count_dir))

        # row index - metadata_avail_cases['BMTcase']
        # BMT_pdMtx = pd.DataFrame(data=bmt_gt_count_dir, index=metadata_avail_cases['BMTcase'], columns = bim['snp'])

        # Remove all duplicated columns (SNPs)
        if len(bim['snp'][bim['snp'].duplicated(False)]) > 0:
            print('>>>> Column(s) {0} have duplicates! Dropping all duplicated columns (SNPs)'.format(list(set(bim['snp'][bim['snp'].duplicated(False)]))))

        BMT_pdMtx = pd.DataFrame(data=BMT_mm_table_count_dir[:, ~(bim['snp'].duplicated(False))],
                                 index=metadata_avail_cases['BMTcase'],
                                 columns=bim['snp'][~(bim['snp'].duplicated(False))])


        BMT_pdMtx.to_hdf(output_fp+'EncodedMatrix/chr'+str(chrom)+'_EncodedMatrix_original_' + mode + '.h5',
                         key='chr_'+str(chrom), complib='blosc', complevel=9)

        # filter 95% call rate
        filtered_index = (1 - BMT_pdMtx.isnull().sum()/BMT_pdMtx.shape[0]) >= 0.95
        BMT_pdMtx_filtered = BMT_pdMtx.loc[:, filtered_index]
        print('Chromosome {0} removed {1} ({2:.2f} %) out of {3} variants (95 % call rate)'.format(chrom,
                                                                                                   (filtered_index == False).sum(),
                                                                                                   (filtered_index == False).sum()/filtered_index.count()*100,
                                                                                                   filtered_index.count()))
        BMT_pdMtx_filtered.to_hdf(output_fp + 'EncodedMatrix/chr' + str(chrom) + '_EncodedMatrix_95filtered_' + mode + '.h5',
                           key='chr_' + str(chrom), complib='blosc', complevel=9)

        # logistic regression
        Encoded_mat = BMT_pdMtx_filtered

        y1 = [metadata_pd.loc[GroupID, 'agvhi24'] for GroupID in Encoded_mat.index]
        y2 = [metadata_pd.loc[GroupID, 'agvhi34'] for GroupID in Encoded_mat.index]

        # number of variants
        num_variants = Encoded_mat.shape[1]

        p_value_table1 = pd.Series(1, index=Encoded_mat.columns)
        p_value_table2 = pd.Series(1, index=Encoded_mat.columns)

        #r_y1 = StrVector(y1)  # aGVHD II ~ IV
        #r_y2 = StrVector(y2)  # aGVHD III ~ IV
        for ind in range(num_variants):
            # x_train = Encoded_mat.iloc[:, ind]
            # y_train = pd.DataFrame(y)
            SNP_mat = pd.DataFrame(Encoded_mat.loc[:, Encoded_mat.columns[ind]])
            SNP_mat['aGVHD24'] = y1
            SNP_mat['aGVHD34'] = y2
            SNP_mat_dropna = SNP_mat.dropna() #drop nan

            r_x = Vector(SNP_mat_dropna.iloc[:, 0])

            r_y1 = StrVector(SNP_mat_dropna.iloc[:, 1])
            r_y2 = StrVector(SNP_mat_dropna.iloc[:, 2])
            try:
                p_value_table1.loc[Encoded_mat.columns[ind]] = list(r_logitRegression(r_x, r_y1))[0]

            except rpy2.rinterface.RRuntimeError as Ee:
                print(Encoded_mat.columns[ind] + ' Cannot detect association (aGVHD24)!')
                p_value_table1.loc[Encoded_mat.columns[ind]] = -1

            try:
                p_value_table2.loc[Encoded_mat.columns[ind]] = list(r_logitRegression(r_x, r_y2))[0]
            except rpy2.rinterface.RRuntimeError as Ee:
                print(Encoded_mat.columns[ind] + ' Cannot detect association! (aGVHD34)')
                p_value_table2.loc[Encoded_mat.columns[ind]] = -1

        p_value_table1.to_hdf(
            output_fp + 'p_values/chr' + str(chrom) + '_logitRegression_p_values_agvhd24_' + mode + '.h5',
            key='chr_' + str(chrom), complib='blosc', complevel=9)

        p_value_table2.to_hdf(
            output_fp + 'p_values/chr' + str(chrom) + '_logitRegression_p_values_agvhd34_' + mode + '.h5',
            key='chr_' + str(chrom), complib='blosc', complevel=9)

        gc.collect()

    except FileNotFoundError as e:
        print(e)







