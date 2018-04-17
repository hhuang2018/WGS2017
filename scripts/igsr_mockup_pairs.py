#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 09:49:28 2018

@author: hhuang2
"""

import pandas as pd
import numpy as np

data_fp = '../Data/'

igsr_tb = pd.read_table(data_fp+'igsr_samples_phase1.tsv')

num_groups = int(igsr_tb.shape[0]/2)

np.random.seed(1)
indices = np.random.permutation(igsr_tb.shape[0])
np.random.seed(1)
CaseIDs = np.random.permutation(range(100000, 1000000))


paired_tb = {'CaseID': ['']*num_groups,
             'DID': ['']*num_groups,
             'RID': ['']*num_groups}
for index in range(num_groups):
    paired_tb['CaseID'][index] = 'IGSR'+str(CaseIDs[index])
    paired_tb['DID'][index] = igsr_tb.loc[indices[index], 'Sample name']
    rid = igsr_tb.loc[indices[igsr_tb.shape[0]-index-1], 'Sample name']
    if rid != paired_tb['DID'][index]:
        paired_tb['RID'][index] = rid
    else:
        print("Redraw the IDs for "+paired_tb['DID'][index])
        
paired_tb_df = pd.DataFrame.from_dict(paired_tb)
paired_tb_df.to_csv(data_fp+'igsr_phase1_PairedIDs.tsv', sep='\t', index = False)

aa = paired_tb['DID'] + paired_tb['RID']
bb = paired_tb_df['DID'].isin(aa)
bbId = bb.index[bb == False].tolist()
cc = paired_tb_df['RID'].isin(aa)
ccId = cc.index[cc == False].tolist()

dd = paired_tb_df['DID'].tolist()
ee = paired_tb_df['RID'].tolist()
ff = dd+ee

aa.sort()
ff.sort()

aa == ff
# len(set(aa))

### HLA typing list
allele_list_fp = '../../2018_CDSW_test/Allelelist.3310.txt'

allele_list = pd.read_table(allele_list_fp, sep = " ") # 17874*2

hla_tb = {'CaseID': ['']*num_groups,
             'DID': ['']*num_groups,
             'RID': ['']*num_groups,
             'r_A_type1_gl': ['']*num_groups,
             'r_A_type2_gl':['']*num_groups,
             'r_B_type1_gl': ['']*num_groups,
             'r_B_type2_gl':['']*num_groups,
             'r_C_type1_gl': ['']*num_groups,
             'r_C_type2_gl':['']*num_groups,
             'r_DRB1_type1_gl': ['']*num_groups,
             'r_DRB1_type2_gl':['']*num_groups,
             'r_DRB3_type1_gl': ['']*num_groups,
             'r_DRB3_type2_gl':['']*num_groups,
             'r_DRB4_type1_gl': ['']*num_groups,
             'r_DRB4_type2_gl':['']*num_groups,
             'r_DRB5_type1_gl': ['']*num_groups,
             'r_DRB5_type2_gl':['']*num_groups,
             'r_DQB1_type1_gl': ['']*num_groups,
             'r_DQB1_type2_gl':['']*num_groups,
             'r_DQA1_type1_gl': ['']*num_groups,
             'r_DQA1_type2_gl':['']*num_groups,
             'r_DPB1_type1_gl': ['']*num_groups,
             'r_DPB1_type2_gl':['']*num_groups,
             'r_DPA1_type1_gl': ['']*num_groups,
             'r_DPA1_type2_gl':['']*num_groups, #### recipient
             'd_A_type1_gl': ['']*num_groups,
             'd_A_type2_gl':['']*num_groups,
             'd_B_type1_gl': ['']*num_groups,
             'd_B_type2_gl':['']*num_groups,
             'd_C_type1_gl': ['']*num_groups,
             'd_C_type2_gl':['']*num_groups,
             'd_DRB1_type1_gl': ['']*num_groups,
             'd_DRB1_type2_gl':['']*num_groups,
             'd_DRB3_type1_gl': ['']*num_groups,
             'd_DRB3_type2_gl':['']*num_groups,
             'd_DRB4_type1_gl': ['']*num_groups,
             'd_DRB4_type2_gl':['']*num_groups,
             'd_DRB5_type1_gl': ['']*num_groups,
             'd_DRB5_type2_gl':['']*num_groups,
             'd_DQB1_type1_gl': ['']*num_groups,
             'd_DQB1_type2_gl':['']*num_groups,
             'd_DQA1_type1_gl': ['']*num_groups,
             'd_DQA1_type2_gl':['']*num_groups,
             'd_DPB1_type1_gl': ['']*num_groups,
             'd_DPB1_type2_gl':['']*num_groups,
             'd_DPA1_type1_gl': ['']*num_groups,
             'd_DPA1_type2_gl':['']*num_groups,
             'agvhi24':['']*num_groups,
             'agvhi34':['']*num_groups}

loci = ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQB1', 'DQA1', 'DPB1', 'DPA1']
#shape, scale = 2. , 2. # gamma distribtuion 
for locus in loci:
        
    locus_allele_list = allele_list[allele_list['AlleleName'].str.contains('^'+locus+'\\*')]
    #denom_sum = locus_allele_list.shape[0] * (locus_allele_list.shape[0]+1)/2
    
    total_num = locus_allele_list.shape[0]
    scale_factor = 10**(len(str(total_num))-1)
    if scale_factor <= 100:
        r = 2
    else:
        r =5
    
    prob = [(total_num+1-pid)*r**(int(total_num/((pid+1)*scale_factor))) for pid in range(total_num)]
    prob = [p/sum(prob) for p in prob]
        
    for index in range(num_groups):
        hla_tb['CaseID'][index] = paired_tb['CaseID'][index]
        hla_tb['DID'][index] = paired_tb['DID'][index]
        hla_tb['RID'][index] = paired_tb['RID'][index]
    
        #np.random.permutation(locus_allele_list.shape[0])
        #np.random.gamma(shape, scale, 1)
        #np.random.exponential(2., 1)
        #aa = np.random.beta(1, 3, 2)
        chosen_ids = np.random.choice(locus_allele_list.shape[0], 2, p=prob)
        chosen_ids.sort()
        
        if np.random.randn() >= 0: # homozygous
            random_assignment = locus_allele_list.iloc[chosen_ids[0]]['AlleleName']
            for ps in ['1', '2']:
                hla_tb['r_'+locus+'_type'+ps+'_gl'][index] = random_assignment
                hla_tb['d_'+locus+'_type'+ps+'_gl'][index] = random_assignment
            
        else: # heterozygous
            for psID in [0, 1]:
                ps = str(psID+1)
                random_assignment = locus_allele_list.iloc[chosen_ids[psID]]['AlleleName']
            
                hla_tb['r_'+locus+'_type'+ps+'_gl'][index] = random_assignment
                hla_tb['d_'+locus+'_type'+ps+'_gl'][index] = random_assignment
        
        if locus == 'A': # outcome
            if np.random.randn() >= 0: # gvhd24
                hla_tb['agvhi24'][index] = 1
                if np.random.randn() >= 0: # gvhd24
                    hla_tb['agvhi34'][index] = 1
                else:
                    hla_tb['agvhi34'][index] = 0
            else:
                hla_tb['agvhi24'][index] = 0
                hla_tb['agvhi34'][index] = 0
                    

hla_tb_df = pd.DataFrame.from_dict(hla_tb)
hla_tb_df.to_csv(data_fp+'igsr_phase1_PairedIDs_HLAtyping-outcome.tsv', sep='\t', index = False)

