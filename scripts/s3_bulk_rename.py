#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
# import glob
import os

s3_bucket = 's3://nmdp-human-longevity-data/modified/renamed_vcf/'
HLI_metadata_fp = '/home/hhuang/data/HLI_metadata/ID_table_wIndex.csv'

# all_vcf_files = os.system('aws s3 ls '+s3_bucket+' --exclude "*" --include  "*.vcf.gz" ')

ID_list = pd.read_csv(HLI_metadata_fp)

num_files = len(ID_list)

for idx in range(num_files):

    original_vcf_prefix = ID_list.Group[idx] + '_' + str(ID_list.GroupID[idx]) + '_' + \
                          ID_list.subjectType[idx] + '_' + str(ID_list.R_D_ID[idx])

    new_vcf_prefix = 'SEQ' + str(ID_list.SeqID[idx])

    os.system('aws s3 mv ' + s3_bucket + original_vcf_prefix + '.vcf.gz ' + s3_bucket + new_vcf_prefix + '.vcf.gz')
    os.system('aws s3 mv ' + s3_bucket + original_vcf_prefix + '.vcf.gz.tbi ' + s3_bucket + new_vcf_prefix + '.vcf.gz.tbi')