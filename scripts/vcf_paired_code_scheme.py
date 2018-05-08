#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 14:50:55 2018

@author: hhuang2
"""


import allel # scikit-allel
import zarr
import numcodecs
import numpy as np
import pandas as pd
import sys
## testing scikit-allel
#allel.__version__
#zarr.__version__
#numcodecs.__version__
#np.__version__


## Data source
ID_table = pd.read_csv('../Data/match_combined_filtered/ID_table.csv')
VCFdata_fp = "../Data/match_combined_filtered/all_chr22.vcf.gz"

zarr_path = '../Data/match_combined_filtered/'

# !ls -lh {VCFdata_fp} # list the files in the directory

# format conversion
allel.vcf_to_zarr(VCFdata_fp, 
                  zarr_path+'all_chr22.zarr', 
                  group='chr22', fields='*', log=sys.stdout, overwrite=True)

callset = zarr.open_group(zarr_path+'all_chr22.zarr', mode='r')
callset.tree(expand=True)

gt_zarr = callset['chr22/calldata/GT']
gt_zarr.info

#pos = callset['chr22/variants/POS']

#loc_region = pos.locate_range(20000, 20100)
samples = callset['chr22/samples'] # columns
gt = allel.GenotypeArray(gt_zarr)  # genotypes

### 
# callset2 = allel.read_vcf(VCFdata_fp)



















