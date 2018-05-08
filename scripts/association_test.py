#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:31:17 2018

@author: hhuang2
"""

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns
#sns.set_style('white')
#sns.set_style('ticks')
#sns.set_context('notebook')
import zarr
#import h5py
import allel; print('scikit-allel', allel.__version__)

from optparse import OptionParser

# plot variant density in windows over the chromosome.
def plot_windowed_variant_density(pos, window_size, title=None):
    # setup windows
    bins = np.arange(0, pos.max(), window_size)

    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1]) / 2

    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size

    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
    return fig


# a function to plot a frequency distribution for any variant attribute.
def plot_variant_hist(f, bins=30):
    x = variants[f][:]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('No. variants')
    ax.set_title('Variant %s distribution' % f)
    plt.gca().set_xticks([2, 3, 4])
    return fig

#  functions to compare two genotypes and return an encoded value
def encode_mmCount(genotype1):  # , genotype2 = None):
    # if genotype2 == None:
    temp = genotype1[:2]
    genotype2 = genotype1[2:]
    genotype1 = temp

    if -1 in genotype1 or -1 in genotype2:  # missing genotypes
        return -1
    elif len(set(genotype1) & set(genotype2)) == 2:  # both matched
        return 0
    elif len(set(genotype1) & set(genotype2)) == 1:  # one matched
        return 1
    elif len(set(genotype1) & set(genotype2)) == 0:  # both mismatched
        return 2


def encode_mmPresAbs(genotype1):  # , genotype2 = None):
    # if genotype2 == None:
    temp = genotype1[:2]
    genotype2 = genotype1[2:]
    genotype1 = temp

    if -1 in genotype1 or -1 in genotype2:  # missing genotypes
        return -1
    elif len(set(genotype1) & set(genotype2)) == 2:  # both matched
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
        #row_id.extend([(pos, 'total')])
    return(row_id)


def contingency_table(encoded_matrix, sample_table):
    '''
    For each position, create a contingency table.
    '''
    caseIDs = sample_table[sample_table.Group == 'a'].GroupID.unique()
    controlIDs = sample_table[sample_table.Group == 'n'].GroupID.unique()

    if encoded_matrix.shape[0] > 1:
        # multi row contingency table
        row_index = add_case_control(
            encoded_matrix.index)  # [(iid, 'case'), (iid, 'control'), (iid, 'total') for iid in row_index]
        caseTable = encoded_matrix.loc[:, [str(iid) for iid in caseIDs]]
        controlTable = encoded_matrix.loc[:, [str(iid) for iid in controlIDs]]
    else:
        # one position table
        row_index = None
        caseTable = encoded_matrix.loc[[str(iid) for iid in caseIDs]]
        controlTable = encoded_matrix.loc[[str(iid) for iid in controlIDs]]

    if row_index != None:
        conting_table = pd.DataFrame({'0': 0,
                                      '1': 0,
                                      '2': 0},
                                     index=pd.MultiIndex.from_tuples(
                                         row_index,
                                         names=['POS', 'case-control']))
        p_value_table = pd.Series(0, index=encoded_matrix.index)

        for pos in encoded_matrix.index:
            caseSumTable = caseTable.loc[pos, :].value_counts()
            controlSumTable = controlTable.loc[pos, :].value_counts()
            for index1 in caseSumTable.index:
                conting_table.loc[pos, 'case'][index1] = caseSumTable[index1]
            for index2 in controlSumTable.index:
                conting_table.loc[pos, 'control'][index2] = controlSumTable[index2]
            # conting_table.loc[pos,'total'] = np.sum(conting_table.loc[pos,])
            try:
                _, p_value_table.loc[pos], _, _ = sp.stats.chi2_contingency(conting_table.loc[pos,])
            except ValueError as ve:
                p_value_table.loc[pos] = -1
    else:
        conting_table = pd.DataFrame({'0': 0,
                                      '1': 0,
                                      '2': 0},
                                     index=['case', 'control'])
        p_value_table = 0

        for index1 in caseSumTable.index:
            conting_table.loc[pos, 'case'][index1] = caseSumTable[index1]
        for index2 in controlSumTable.index:
            conting_table.loc[pos, 'control'][index2] = controlSumTable[index2]

        try:
            _, p_value_table, _, _ = sp.stats.chi2_contingency(conting_table)
        except ValueError as ve:
            p_value_table = -1

    return conting_table, p_value_table

# plot p-values
def plot_p_value(p_value_table, title, y_label):
    '''
    plot p-values if it's an array
    '''
    if len(p_value_table) > 1:
        fig, ax = plt.subplots(figsize=(12, 4))
        sns.despine(ax=ax, offset=10)
        left = np.arange(len(p_value_table))
        palette = sns.color_palette()
        pop2color = {'a': palette[0], 'n': palette[0]}
        colors = [pop2color[p] for p in samples_subset.Group]
        ax.bar(left, p_value_table, color=colors)
        ax.set_xlim(0, len(p_value_table))
        ax.set_xlabel('Posisitons')
        ax.set_ylabel(y_label)
        ax.set_title(title)
        plt.axhline(y=yline, linestyle="dashdot", color='r')
        #handles = [mpl.patches.Patch(color=palette[0]),
        #           mpl.patches.Patch(color=palette[1])]
        return fig

### File directory
#samples_fn = '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/ID_table_wIndex.csv'

#output_fn = '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/'

#zarr_path = '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/all_by_chr/'
#callset = zarr.open_group(zarr_path+'all_chr22_genotypes.zarr', mode='r')
# callset.tree(expand=True)


### parse input
parser = OptionParser()
#parser.add_option("-i", "--input", dest="input_filename",
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
                  help="Chromosome ID")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

zarr_path = options.input_filename
chromID = options.chromID
samples_fn = options.sample_metadata
output_fp = options.output_fp

## load variants
callset = zarr.open_group(zarr_path, mode='r')
variants = allel.VariantChunkedTable(callset['chr'+chromID]['variants'],
                                     names=['POS', 'REF', 'ALT', 'AN', 'AC', 'numalt'],
                                     index='POS')
#chrom = 'chr22'
#variants = allel.VariantChunkedTable(callset[chrom]['variants'],
#                                     names=['POS', 'REF', 'ALT', 'AN', 'AC', 'numalt'],
#                                     index='POS')
# variants

# POS
pos = variants['POS'][:]
# pos

# a plot with the SNP positions from our chosen chromosome.
f = plot_windowed_variant_density(pos, window_size=100000, title='Raw variant density')
f.savefig(output_fp+'/'+chromID+'RawVariantsDensity.pdf', bbox_inches='tight')

# See how many biallelic, triallelic and quadriallelic SNPs we have.
f2 = plot_variant_hist('numalt', bins=np.arange(1.5, 5.5, 1))
# plt.gca().set_xticks([2, 3, 4])
f2.savefig(output_fp+'/'+chromID+'MultiAlelles.pdf', bbox_inches='tight')

#### variant filtering
filter_expression = '(AN >= 800)'
variant_selection = variants.eval(filter_expression)[:]
# variant_selection

print("Kept {} variants and removed {} variants".format(np.count_nonzero(variant_selection), np.count_nonzero(~variant_selection)))

variants_pass_pos = pos.compress(variant_selection)
# variants_pass_pos

## Genoyptes
genotypes = allel.GenotypeChunkedArray(callset['chr'+chromID+'/calldata/GT'])
genotypes_subset = genotypes.subset(variant_selection, )

#### Sample metadata
# samples_fn = '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/ID_table_wIndex.csv'
samples = pd.DataFrame.from_csv(samples_fn)
# samples.head()

# paired groups filtering
paired_groupIDs_list = samples.GroupID.value_counts() == 2 # D-R pairs
#aa = paired_groupIDs.to_frame()
paired_groupIDs = paired_groupIDs_list.loc[paired_groupIDs_list==True].index
# paired_groupIDs
sample_selection = samples.GroupID.isin(paired_groupIDs).values
samples_subset = samples[sample_selection]
samples_subset.reset_index(drop=True, inplace=True)
# samples_subset.head()

Unique_GroupIDs = samples_subset.GroupID.unique()

##### Encoding - scheme 1: count-based
encodedMatrix = pandas.DataFrame(data=-1, columns=[str(ggid) for ggid in Unique_GroupIDs],
                                 index=[str(posID) for posID in variants_pass_pos])
cc = 0
for gpid in Unique_GroupIDs:
    print("counter = {}; gpid= {}".format(cc, gpid))
    # gpid = Unique_GroupIDs[ind]
    IDs = samples_subset[samples_subset.GroupID == gpid]
    donor_index = SampleIDs.index(str(IDs[IDs.subjectType == 'D'].SeqID.item()))
    recipient_index = SampleIDs.index(str(IDs[IDs.subjectType == 'R'].SeqID.item()))
    combined_genotypes = np.hstack((genotypes_subset[:, donor_index, :], genotypes_subset[:, recipient_index, :]))

    encodedMatrix.loc[:, str(gpid)] = np.apply_along_axis(encode_mmCount, 1, combined_genotypes)
    cc += 1

###### Encoding - scheme 2: presence/absence-based
encodedMatrix_presAbs = pandas.DataFrame(data=-1, columns=[str(ggid) for ggid in Unique_GroupIDs],
                                         index=[str(posID) for posID in variants_pass_pos])
cc = 0
for gpid in Unique_GroupIDs:
    print("counter = {}; gpid= {}".format(cc, gpid))
    # gpid = Unique_GroupIDs[ind]
    IDs = samples_subset[samples_subset.GroupID == gpid]
    donor_index = SampleIDs.index(str(IDs[IDs.subjectType == 'D'].SeqID.item()))
    recipient_index = SampleIDs.index(str(IDs[IDs.subjectType == 'R'].SeqID.item()))
    combined_genotypes = np.hstack((genotypes_subset[:, donor_index, :], genotypes_subset[:, recipient_index, :]))

    encodedMatrix_presAbs.loc[:, str(gpid)] = np.apply_along_axis(encode_mmPresAbs, 1, combined_genotypes)
    cc += 1

##### chi-square test - count-based scheme
alpha = 0.05
contin_table, p_value_tb_count = contingency_table(encodedMatrix, samples_subset)
f3 = plot_p_value(p_value_tb_count, 'Count-based Scheme Chi-square test p-values', alpha, 'p-values')
f3.savefig(output_fp+'/'+chromID+'CountBased_pvalues.pdf', bbox_inches='tight')

positive_p_values_count = p_value_tb_count.loc[p_value_tb_count > 0]
negative_p_values_count = p_value_tb_count.loc[p_value_tb_count <= 0]

log_p_values_count = -np.log(positive_p_values_count)
f4 = plot_p_value(log_p_values_count, 'Count-based Scheme Chi-square test p-values (-log(p))', -np.log(alpha), '-log(p)')
f4.savefig(output_fp+'/'+chromID+'CountBased_Log_pvalues.pdf', bbox_inches='tight')

#### chi-square test - presence/absence-based scheme
contin_table, p_value_tb_presabs = contingency_table(encodedMatrix, samples_subset)
f5 = plot_p_value(p_value_tb_presabs, 'Presence/absence-based Scheme Chi-square test p-values', alpha, 'p-values')
f5.savefig(output_fp+'/'+chromID+'Presabsence_pvalues.pdf', bbox_inches='tight')

positive_p_values_presabs = p_value_tb_presabs.loc[p_value_tb_presabs > 0]
negative_p_values_presabs = p_value_tb_presabs.loc[p_value_tb_presabs <= 0]

log_p_values_presabs = -np.log(positive_p_values_presabs)
f6 = plot_p_value(log_p_values_presabs, 'Presence/absence-based Scheme Chi-square test p-values (-log(p))', -np.log(alpha), '-log(p)')
f6.savefig(output_fp+'/'+chromID+'Presabsence_Log_pvalues.pdf', bbox_inches='tight')