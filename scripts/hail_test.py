#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 09:21:11 2018

@author: hhuang2
"""

import hail
hc = hail.HailContext()

for chrom in range(1, 23):
    hc.import_plink(bed = 's3://nmdp-plink-bucket/plink-files/ImputeQC/'+ str(chrom) +'bed.bed',
                    bim = 's3://nmdp-plink-bucket/plink-files/ImputeQC/'+ str(chrom) +'bed.bim',
                    fam = 's3://nmdp-plink-bucket/plink-files/ImputeQC/'+ str(chrom) +'bed.fam').write('nmdp-hhuang/GWAS_dev/Hail_VDS/hail_chr'+chrom+'.vds')
    
