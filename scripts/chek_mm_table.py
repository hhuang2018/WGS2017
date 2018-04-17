#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 19:16:55 2018

@author: hhuang2
"""

import sys
sys.path.append('../utils/')
from utils import coreFunctions as cf
from pandas_plink import read_plink
#from pandas import DataFrame
import numpy as np
import feather

metadata_fp = "../Data/available_cases.csv"
output = '../Data/'

chrom = 22
BMT_mm_table_n = np.load(output+'BMT_mm_table_chr_'+str(chrom)+'.npz')

mm_table = BMT_mm_table_n['mm_table']
ID_table = BMT_mm_table_n['ID_table']

feather.write_dataframe(mm_table, output+'mm_table12_test.feather')
