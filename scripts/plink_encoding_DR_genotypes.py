#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur July 19 00:00:00 2018

@author: hhuang2
"""
from utils import coreFunctions as cf
import numpy as np
import pandas as pd

from pandas_plink import read_plink

import gc

mode = 'DR_genotype_encoding'
## EC2
plink_fp = "/home/hhuang/efs/GWASH_IMPUTED_DATA/ImputeQC/"
metadata_fp = "/home/hhuang/data/GWAS/available_cases.csv"
output_fp = '/home/hhuang/data/GWAS/MismatchEncoded/' + mode + '/'

chrList = range(1, 23)  # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])

metadata_pd = pd.read_csv(metadata_fp, index_col=2)
metadata_pd.index = metadata_pd.index.map(str)

for chrom in chrList:
    # BMT_mm_table_count_dir = np.array([], dtype='float64')

    try:
        (bim, fam, bed) = read_plink(plink_fp + str(chrom) + 'bed')
        # bim - pandas.DataFrame – Alleles.
        # fam - pandas.DataFrame – Samples.
        # G - numpy.ndarray – Genotype.

        DR_gt_encoding = np.empty([num_case, bim.shape[0]], dtype='float64') - 1

        # fam.iid.head()
        # G.head()

        for case_index in range(num_case):

            bmt_fams = fam.query("iid in ['" + metadata_avail_cases['NMDP_DID'][case_index] +
                                 "', '" +
                                 metadata_avail_cases['NMDP_RID'][case_index] + "']")
            gt = bed[:, bmt_fams.i.values].compute()

            bmt_gt = gt[:, 0] * 3 + gt[:, 1]  # donor*3+recipient = gt code

            DR_gt_encoding[case_index, :] = bmt_gt

        # Remove all duplicated columns (SNPs)
        if len(bim['snp'][bim['snp'].duplicated(False)]) > 0:
            print('>>>> Column(s) {0} have duplicates! Dropping all duplicated columns (SNPs)'.format(list(set(bim['snp'][bim['snp'].duplicated(False)]))))

        BMT_pdMtx = pd.DataFrame(data=DR_gt_encoding[:, ~(bim['snp'].duplicated(False))],
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

