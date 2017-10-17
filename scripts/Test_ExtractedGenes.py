#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 22:16:54 2017

@author: hhuang2
"""

import glob
import csv

file_fp = "../HLI/GeneData/"

GeneName = "IL10"
Gene_files = glob.glob(file_fp+"*_"+GeneName+".vcf")


CHR=[]
POS=[]
#f = open(Gene_files[0],'r') 
with open(Gene_files[0], 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        print(row)