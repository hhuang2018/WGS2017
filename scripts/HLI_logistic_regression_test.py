# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:31:17 2018

@author: hhuang2
"""
import pickle
import numpy as np
import scipy as sp
import pandas as pd
import rpy2
from rpy2.robjects.packages import importr
#from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro

#import zarr

from optparse import OptionParser


# R functions
stats = importr('stats')
base = importr('base')
#Matrix = ro.r.matrix
Vector = rpy2.robjects.vectors.IntVector
StrVector = rpy2.robjects.vectors.StrVector
ro.r('''
        # create a function `LogisticRegression`
        LogitRegression <- function(x, y) {
            x[x==-1] = NA
            snpData = cbind(x, y)
            
            snpData = as.data.frame(snpData)
            colnames(snpData) = c('gt', 'Outcome')
            
            Logit_regression <- glm(formula = Outcome ~ gt, family = binomial(), data = snpData)
            
            p_val <- coef(summary(Logit_regression))[2,'Pr(>|z|)']
            return(p_val)
            # return(summary(regression))
        }
        # call the function `f` with argument value 3
        
        ''')

r_logitRegression = ro.globalenv['LogitRegression']

### parse input
parser = OptionParser()

parser.add_option("-s", "--sample", dest="sample_metadata",
                  action='store', type="string",
                  help="Sample metatdata", metavar="FILE")
parser.add_option("-e", "--encodedmatrix", dest="encoded_mat_filename",
                  action='store', type="string",
                  help="Encoded matrix file", metavar="FILE")
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

encoded_mat_filename = options.encoded_mat_filename
chromID = options.chromID
samples_fn = options.sample_metadata
output_fp = options.output_fp

## load encoded matrix

with open(encoded_mat_filename, 'rb') as eM:
    EncodedMatrix = pickle.load(eM)

# Load sample metadata
samples = pd.read_csv(samples_fn)
labels = samples[['Group', 'GroupID']]

# encoded matrix - sampleID * Location
Encoded_mat = EncodedMatrix.T

# label
y = [list(labels.Group[labels.GroupID == int(GroupID)])[0] for GroupID in Encoded_mat.index]


# number of variants
num_variants = Encoded_mat.shape[1]

p_value_table = pd.Series(1, index=Encoded_mat.columns)

# test each varaint
r_y = StrVector(y)
for ind in range(num_variants):
    #x_train = Encoded_mat.iloc[:, ind]
    #y_train = pd.DataFrame(y)

    r_x = Vector(Encoded_mat.loc[:, Encoded_mat.columns[ind]])
    try:
        p_value_table.loc[Encoded_mat.columns[ind]] = list(r_logitRegression(r_x, r_y))[0]
    except rpy2.rinterface.RRuntimeError as Ee:
        print(Encoded_mat.columns[ind] + ' no association!')

p_value_table.to_hdf(output_fp+'chr'+chromID+'_logitRegression_p_values.h5', key='chr_'+chromID, complib='blosc', complevel=9)

