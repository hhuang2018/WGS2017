#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:31:17 2018

@author: hhuang2
"""
import glob
import numpy as np
from scipy import stats as sci_stats
import scipy as sp
import pandas as pd
import re
import pickle

#  functions to compare two genotypes and return an encoded value
def encode_mmCount(genotype1):  # , genotype2 = None):
    # if genotype2 == None:
    temp = genotype1[:2]
    genotype2 = genotype1[2:]
    genotype1 = temp

    common_element_len = len(set(genotype1) & set(genotype2))
    new_code = -1
    if -1 in genotype1 or -1 in genotype2:  # missing genotypes
        new_code = -1
    elif common_element_len == 2:  # len(set(genotype1) & set(genotype2)) == 2:  # both matched
        new_code = 0
    elif common_element_len == 1:  # one matched
        if len(set(genotype1)) == 1 and len(set(genotype2)) == 1:
            new_code = 0  # both matched
        else:
            new_code = 1  # one mismatch
    elif common_element_len == 0:  # both mismatched
        new_code = 2

    return new_code


def encode_mmPresAbs(genotype1):  # , genotype2 = None):
    # if genotype2 == None:
    temp = genotype1[:2]
    genotype2 = genotype1[2:]
    genotype1 = temp

    if -1 in genotype1 or -1 in genotype2:  # missing genotypes
        return -1
    elif len(set(genotype1) & set(genotype2)) == 2:  # both matched
        return 0
    elif len(set(genotype1) & set(genotype2)) == 1 and len(set(genotype1)) == 1 and len(set(genotype2)) == 1:
        return 0
    else:  # one or two mismatched
        return 1


# chi-square test and p-values
def add_case_control(pos_list):
    '''
    Helper function to create multi-layered row index
    '''
    row_id = []
    for pos in pos_list:
        row_id.extend([(pos, 'case')])
        row_id.extend([(pos, 'control')])
        # row_id.extend([(pos, 'total')])
    return (row_id)


def contingency_table1(encoded_matrix, sample_table, mode='count'):
    '''
    For each position, create a contingency table.
    '''
    caseIDs = sample_table[sample_table.Group == 'a'].GroupID.unique()
    controlIDs = sample_table[sample_table.Group == 'n'].GroupID.unique()
    uniq_index = np.unique(encoded_matrix.index)

    if encoded_matrix.shape[0] > 1:
        # multi row contingency table

        row_index = add_case_control(
            uniq_index)  # [(iid, 'case'), (iid, 'control'), (iid, 'total') for iid in row_index]
        caseTable = encoded_matrix.loc[:, [str(iid) for iid in caseIDs]]
        controlTable = encoded_matrix.loc[:, [str(iid) for iid in controlIDs]]
    else:
        # one position table
        row_index = None
        caseTable = encoded_matrix.loc[[str(iid) for iid in caseIDs]]
        controlTable = encoded_matrix.loc[[str(iid) for iid in controlIDs]]

    if row_index != None:
        if mode[0] == 'c':
            conting_table = pd.DataFrame({'0': 1,
                                          '1': 1,
                                          '2': 1},
                                         index=pd.MultiIndex.from_tuples(
                                             row_index,
                                             names=['POS', 'case-control']))
        elif mode[0] == 'p':
            conting_table = pd.DataFrame({'0': 1,
                                          '1': 1},
                                         index=pd.MultiIndex.from_tuples(
                                             row_index,
                                             names=['POS', 'case-control']))
        p_value_table = pd.Series(0, index=uniq_index)

        for pos in uniq_index:

            try:
                caseSumTable = caseTable.loc[pos, :].value_counts()
            except AttributeError as ae:  # duplicated location
                temp_tb = caseTable.loc[pos, :].sum()
                caseSumTable = temp_tb.astype('int').value_counts()
            try:
                controlSumTable = controlTable.loc[pos, :].value_counts()
            except AttributeError as ae:
                temp_tb = controlTable.loc[pos, :].sum()
                controlSumTable = temp_tb.astype('int').value_counts()

            for index1 in caseSumTable.index:
                if index1 != -1:
                    conting_table.loc[pos, 'case'][index1] = caseSumTable[index1]

            for index2 in controlSumTable.index:
                if index2 != -1:
                    conting_table.loc[pos, 'control'][index2] = controlSumTable[index2]
            # conting_table.loc[pos,'total'] = np.sum(conting_table.loc[pos,])
            try:
                _, p_value_table.loc[pos], _, _ = sp.stats.chi2_contingency(conting_table.loc[pos,])
            except ValueError as ve:
                p_value_table.loc[pos] = -1

    else:
        if mode[0] == 'c':
            conting_table = pd.DataFrame({'0': 1,
                                          '1': 1,
                                          '2': 1},
                                         index=['case', 'control'])
        elif mode[0] == 'p':
            conting_table = pd.DataFrame({'0': 1,
                                          '1': 1},
                                         index=['case', 'control'])
        p_value_table = 0

        for index1 in caseSumTable.index:
            if index1 != -1:
                conting_table.loc[pos, 'case'][index1] = caseSumTable[index1]
        for index2 in controlSumTable.index:
            if index2 != -1:
                conting_table.loc[pos, 'control'][index2] = controlSumTable[index2]

        try:
            _, p_value_table, _, _ = sp.stats.chi2_contingency(conting_table)
        except ValueError as ve:
            p_value_table = -1

    return conting_table, p_value_table

# Load sample  metadata
# samples_fn = '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/ID_table_wIndex.csv'
samples = pd.read_csv('/home/hhuang/data/HLI_metadata/ID_table_wIndex.csv')
# samples.head()

# paired groups filtering
paired_groupIDs_list = samples.GroupID.value_counts() == 2  # D-R pairs
# aa = paired_groupIDs.to_frame()
paired_groupIDs = paired_groupIDs_list.loc[paired_groupIDs_list == True].index
# paired_groupIDs
sample_selection = samples.GroupID.isin(paired_groupIDs).values
samples_subset = samples[sample_selection]
samples_subset.reset_index(drop=True, inplace=True)

# Load encoded matrices - count based
count_based_encodedMatrix_file = glob.glob('chr*EncodedMatrix_CountBased.pkl')

counter = 0
for EnMat_fp in count_based_encodedMatrix_file:

    chromID = re.sub('chr', '', EnMat_fp.split('_')[0])

    with open(EnMat_fp, 'rb') as eM:
        EncodedMatrix = pickle.load(eM)

    # chi-square test - count-based scheme
    # alpha = 0.05
    conting_table_count, p_value_tb_count = contingency_table1(EncodedMatrix, samples_subset, 'count')

    conting_table_count.to_hdf('Count_based_pseudoCount/Contingency_tb_count_' + EnMat_fp.split('_')[0] + '.h5', key=EnMat_fp.split('_')[0])

    if counter == 0:
        all_p_values_count = pd.DataFrame({'POS': p_value_tb_count.index,
                                           'p_value': p_value_tb_count,
                                           'chromosome': EnMat_fp.split('_')[0]})
    else:
        temp_df = pd.DataFrame({'POS': p_value_tb_count.index,
                                'p_value': p_value_tb_count,
                                'chromosome': EnMat_fp.split('_')[0]})

        all_p_values_count = pd.concat([all_p_values_count, temp_df])

    counter += 1


all_p_values_count.to_hdf('Count_based_pseudoCount/All_chr_p_values_countbased.h5', key= 'all_chr_p', complib='blosc', complevel=9)

all_p_values_count_postive = all_p_values_count.loc[all_p_values_count.p_value >= 0]
all_p_values_count_postive.to_csv('Count_based_pseudoCount/All_chr_positive_pvals_countbased.csv')


# load pres/abs-based encoded matrix
count_based_encodedMatrix_file = glob.glob('chr*EncodedMatrix_presabs.pkl')

counter = 0
for EnMatPres_fp in count_based_encodedMatrix_file:

    #chromID = re.sub('chr', '', EnMatPres_fp.split('_')[0])

    with open(EnMatPres_fp, 'rb') as eM_p:
        EncodedMatrix_pres = pickle.load(eM_p)

    # chi-square test - count-based scheme
    # alpha = 0.05
    conting_table_presabs, p_value_tb_presabs = contingency_table1(EncodedMatrix_pres, samples_subset, 'presence')

    conting_table_presabs.to_hdf('PresAbs_based_pseudoCount/Contingency_tb_PresAbs_' + EnMatPres_fp.split('_')[0] + '.h5',
                               key=EnMatPres_fp.split('_')[0])

    if counter == 0:
        all_p_values_presabs = pd.DataFrame({'POS': p_value_tb_presabs.index,
                                           'p_value': p_value_tb_presabs,
                                           'chromosome': EnMatPres_fp.split('_')[0]})
    else:
        temp_df = pd.DataFrame({'POS': p_value_tb_presabs.index,
                                'p_value': p_value_tb_presabs,
                                'chromosome': EnMatPres_fp.split('_')[0]})

        all_p_values_presabs = pd.concat([all_p_values_presabs, temp_df])

    counter += 1

all_p_values_presabs.to_hdf('PresAbs_based_pseudoCount/All_chr_p_values_presabs.h5', key='all_chr_p', complib='blosc',
                    complevel=9)

all_p_values_presabs_postive = all_p_values_presabs.loc[all_p_values_presabs.p_value >= 0]
all_p_values_presabs_postive.to_csv('PresAbs_based_pseudoCount/All_chr_positive_pvals_PresAbs.csv')

