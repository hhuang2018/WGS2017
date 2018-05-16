#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:31:17 2018

@author: hhuang2
"""

import numpy as np
from scipy import stats as sci_stats
import scipy as sp
import pandas as pd

import zarr
import allel

from optparse import OptionParser


# Group encoding
def encode_group_mmtype(genotypes):
    genotype1 = genotypes[:2]
    genotype2 = genotypes[2:]

    common_gt = set(genotype1) & set(genotype2)

    if -1 in genotype1 or -1 in genotype2:  # missing genotypes
        new_code = [-1, -1]
    elif len(common_gt) == 2:  # both matched
        new_code = [0, 0]
    elif len(common_gt) == 1:  # one matched
        if len(set(genotype1)) == 1 and len(set(genotype2)) == 1:
            new_code = [0, 0]
        else:
            new_gt1 = remaining_gt(genotype1, common_gt)
            new_gt2 = remaining_gt(genotype2, common_gt)
            new_code = [0, encoded_mmtype(new_gt1, new_gt2)]
    elif len(set(genotype1) & set(genotype2)) == 0:  # both mismatched
        new_code = []
        for ind in range(2):
            new_code.append(encoded_mmtype(genotype1[ind], genotype2[ind]))
    else:
        new_code = [-1, -1]

    return new_code


def remaining_gt(genotype, common_gt):
    if len(set(genotype) - common_gt) == 1:  # mismatched type
        gt = list(set(genotype) - common_gt)[0]
    else:
        gt = list(common_gt)[0]

    return gt


def encoded_mmtype(gt1, gt2):
    '''
    gt1 - one allele from donor
    gt2 - one allele from recipient
    encoding scheme:
    gt1 --> gt2   ---- new_code
     0      1             1
     0      2             2
     1      0             3
     1      2             4
     2      0             5
     2      1             6
    '''

    new_code = -1

    if gt1 == 0:
        if gt2 == 1:
            new_code = 1
        else:
            new_code = 2
    elif gt1 == 1:
        if gt2 == 0:
            new_code = 3
        else:
            new_code = 4
    else:
        if gt2 == 0:
            new_code = 5
        else:
            new_code = 6

    return new_code


def frequency_table(gt_row):
    gt_row = pd.Series(gt_row)

    freq_table = np.zeros(8)
    stats_tb = gt_row.value_counts()

    for i in range(len(stats_tb)):
        freq_table[stats_tb.index[i]] = stats_tb.iloc[i]

    return freq_table

# parse input
parser = OptionParser()

# p arser.add_option("-i", "--input", dest="input_filename",
#                  action = 'store', type = "string",
#                  help="Open Zarr file", metavar="FILE")

parser.add_option("-s", "--sample", dest="sample_metadata",
                  action = 'store', type = "string",
                  help="Sample metatdata", metavar="FILE")
parser.add_option("-z", "--zarrfile", dest="zarr_filename",
                  action = 'store', type = "string",
                  help="Zarrr file", metavar="FILE")
parser.add_option("-g", "--group", dest="chromID",
                  action='store', type="string",
                  help="Chromosome ID")
parser.add_option("-o", "--output", dest="output_fp",
                  action='store', type="string",
                  help="output file path")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

zarr_path = options.zarr_filename
chromID = options.chromID
samples_fn = options.sample_metadata
output_fp = options.output_fp

# load variants
callset = zarr.open_group(zarr_path, mode='r')
variants = allel.VariantChunkedTable(callset[chromID]['variants'],
                                     names=['POS', 'REF', 'ALT', 'AN', 'AC', 'numalt'],
                                     index='POS')

filter_expression = '(AN >= 800)'
variant_selection = variants.eval(filter_expression)[:]

pos = variants['POS'][:]
variants_pass_pos = pos.compress(variant_selection)

genotypes = allel.GenotypeChunkedArray(callset[chromID+'/calldata/GT'])
genotypes_subset = genotypes.subset(variant_selection, )

samples = pd.read_csv(samples_fn)

# paired groups filtering
paired_groupIDs_list = samples.GroupID.value_counts() == 2  # D-R pairs
# aa = paired_groupIDs.to_frame()

paired_groupIDs = paired_groupIDs_list.loc[paired_groupIDs_list==True].index
# paired_groupIDs
sample_selection = samples.GroupID.isin(paired_groupIDs).values
samples_subset = samples[sample_selection]
samples_subset.reset_index(drop=True, inplace=True)
# samples_subset.head()

Unique_GroupIDs = samples_subset.GroupID.unique()

SampleIDs = list(callset[chromID+'/samples'])

# All group encoding
groups = [str(ggid) for ggid in Unique_GroupIDs]
groups = [["{}.1".format(g_), "{}.2".format(g_)] for g_ in groups]
col_name = np.array(groups).reshape(len(groups)*2, )

encodedMatrix = pd.DataFrame(data=-1, columns=col_name, index=[str(posID) for posID in variants_pass_pos])

cc = 0
for gpid in Unique_GroupIDs:
    print("counter = {}; gpid= {}".format(cc, gpid))
    # gpid = Unique_GroupIDs[ind]
    IDs = samples_subset[samples_subset.GroupID == gpid]
    donor_index = SampleIDs.index(str(IDs[IDs.subjectType == 'D'].SeqID.item()))
    recipient_index = SampleIDs.index(str(IDs[IDs.subjectType == 'R'].SeqID.item()))
    combined_genotypes = np.hstack((genotypes_subset[:, donor_index, :], genotypes_subset[:, recipient_index, :]))

    encodedMatrix.loc[:, [str(gpid) + '.1', str(gpid) + '.2']] = np.apply_along_axis(encode_group_mmtype, 1,
                                                                                     combined_genotypes)

    cc += 1

encodedMatrix.to_hdf(output_fp+'/chr'+chromID+'_mismatch_type_encoding.h5', key='chr'+chromID)

# calculate the frequencies

freq_col_name = ['Match', 'REF.ALT1', 'REF.ALT2', 'ALT1.REF', 'ALT1.ALT2', 'ALT2.REF', 'ALT2.ALT1', 'Missing']

all_freq_table1 = np.apply_along_axis(frequency_table, 1, encodedMatrix)

all_freq_table = pd.DataFrame(all_freq_table1, columns=freq_col_name, index=[str(posID) for posID in variants_pass_pos])

all_freq_table.to_hdf(output_fp+'/chr'+chromID+'_mismatch_type_frequency_table.h5', key='chr'+chromID)