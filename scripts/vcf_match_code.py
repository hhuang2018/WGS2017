#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:31:17 2018

@author: hhuang2
"""

import vcf # PyVCF - A Variant Call Format Parser for Python pacakge
# http://pyvcf.readthedocs.io/en/latest/

import allel # scikit-allel
import zarr
import numcodecs
import numpy as np
import sys
# http://scikit-allel.readthedocs.io/en/latest/
# http://alimanfoo.github.io/2018/04/09/selecting-variants.html


VCFdata_fp = "../Data/match_combined_filtered/"

## testing PyVCF package

vcf_reader = vcf.Reader(filename=VCFdata_fp+'a_106_filtered.vcf.gz')
vcf_writer = vcf.Writer(open('/dev/null', 'w'), vcf_reader)
for record in vcf_reader:
    
    original_gt_case = [record.samples[0].gt_type, record.samples[1].gt_type]
    # hom_ref = 0 het = 1 hom_alt = 2
    try:
        np.abs(np.diff(original_gt_case))
    except TypeError as te:
        pass
    
    
    # How does it work with tri-allelic variants?
    




############################
## testing scikit-allel
allel.__version__
zarr.__version__
numcodecs.__version__
np.__version__

# !ls -lh {VCFdata_fp} # list the files in the directory

# format conversion
zarr_path = '../Data/match_combined_filtered/'
allel.vcf_to_zarr(VCFdata_fp+'a_106_filtered.vcf.gz', 
                  zarr_path+'a_106_filtered_genotypes.zarr', 
                  group='106', fields='*', log=sys.stdout, overwrite=True)

callset = zarr.open_group(zarr_path+'a_106_filtered_genotypes.zarr', mode='r')
callset.tree(expand=True)

gt_zarr = callset['106/calldata/GT']
gt_zarr.info

pos = callset['106/variants/POS']

loc_region = pos.locate_range(20000000, 20100000)

gt_region = allel.GenotypeArray(gt_zarr)






