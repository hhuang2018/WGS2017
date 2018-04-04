#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 09:33:21 2018

@author: hhuang2
"""
import csv

def isDonor(NMDP_ID):
    """
    Check if the NMDP_ID is a donor or a recipient
    Donor: 4-4-1
    Recipient: 3-3-1
    """
    if NMDP_ID.index("-") == 4:
        return("D")
    else:
        return("R")

def int2DRid(NMDP_ID, DRtype = "donor"):
    """
    Check if the NMDP_ID is a donor or a recipient
    Donor: 4-4-1
    Recipient: 3-3-1
    """
    if DRtype.upper()[0] == 'D':
        temp_NMDP_ID = '0'*(9-len(str(NMDP_ID))) + str(NMDP_ID)
        new_NMDP_ID = temp_NMDP_ID[:4]+ '-'+ temp_NMDP_ID[4:8] + '-' + temp_NMDP_ID[8:]
    else:
        temp_NMDP_ID = '0'*(7-len(str(NMDP_ID))) + str(NMDP_ID)
        new_NMDP_ID = temp_NMDP_ID[:3]+ '-'+ temp_NMDP_ID[3:6] + '-' + temp_NMDP_ID[6:]
    
    return(new_NMDP_ID)

def DRid2int(NMDP_ID):
    """
    Check if the NMDP_ID is a donor or a recipient
    Donor: 4-4-1
    Recipient: 3-3-1
    """
    
    return([int(NMDP_ID.replace("-","")), isDonor(NMDP_ID)])
    
    
def readCaseInfo(fp, header = True):
    """
    Read BMT case ID into a list
    header - flag true if the table has a header; then will remove the header
    """
    # fp = "../../rawData/SG39_caseID.csv"
    
    caseIDs = {"BMTcase":[], "NMDP_DID":[], "NMDP_RID":[]}
    with open(fp, 'r') as f:
        if header:
            h = next(f) # skip header
            headers = h.split(",")
        reader=csv.reader(f)
        for line in reader:
           if line[headers.index('bmt_case_num')] not in caseIDs["BMTcase"]: # new record
               caseIDs["BMTcase"].append(line[headers.index('bmt_case_num')])
               caseIDs["NMDP_RID"].append(int2DRid(line[headers.index('nmdp_rid')], 'R'))
               caseIDs["NMDP_DID"].append(int2DRid(line[headers.index('nmdp_did')], 'D'))

    return(caseIDs)

def twoFieldHLA(tplist):
    '''
    '''
    new_tplist = []
    for tp in tplist:
        new_tplist.append(":".join(tp.split(":")[:2]))
    new_tplist.sort()
    return(new_tplist)

def readHLAtable(fp, header = True):
    """
    Read BMT case HLA into a list
    header - flag true if the table has a header; then will remove the header
    """
    caseIDs = {"BMTcase":[], "NMDP_DID":[], "NMDP_RID":[], 
               "HLA-A":[], "HLA-B": [], "HLA-C":[],
               "HLA-DRB1":[], "HLA-DQB1":[], "HLA-DPB1":[],
               "transplantDate":[], 
               "Recipient_BroadRace":[], "Recipient_rollupRace":[],
               "Donor_BroadRace":[], "Donor_rollupRace":[]}
    with open(fp, 'r') as f:
        if header:
            h = next(f) # skip header
            headers = h.split(",")
            
        reader=csv.reader(f)
        for line in reader:
           if line[headers.index('bmt_case_num')] not in caseIDs["BMTcase"]: # new record
               caseIDs["BMTcase"].append(line[headers.index('bmt_case_num')])
               caseIDs["NMDP_RID"].append(int2DRid(line[headers.index('nmdp_rid')], 'R'))
               caseIDs["NMDP_DID"].append(int2DRid(line[headers.index('nmdp_did')], 'D'))
               
               # txDate: 85
               caseIDs["transplantDate"].append(line[headers.index('transplant_date')])
               # Recipient BroadRace: 88
               caseIDs["Recipient_BroadRace"].append(line[headers.index('ridbroadrace')])
               # Recipient rollUpRace: 89
               caseIDs["Recipient_rollupRace"].append(line[headers.index('ridrolluprace')])
               # Donor BroadRace: 91
               caseIDs["Donor_BroadRace"].append(line[headers.index('didbroadrace')])
               # Donor rollUpRace: 92
               caseIDs["Donor_rollupRace"].append(line[headers.index('didrolluprace')])
               
               # HLA-A: 4,5, 26, 27
               r_type = [line[headers.index('r_a_typ1')], line[headers.index('r_a_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_a_typ1')], line[headers.index('d_a_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-A"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-A"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-A"].append('MM_'+'_'.join(r_type))
                   #caseIDs["HLA-A"].append('MM_'+'_'.join(r_type.extend(d_type)))
                   
               # HLA-B: 6,7, 28, 29
               r_type = [line[headers.index('r_b_typ1')], line[headers.index('r_b_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_b_typ1')], line[headers.index('d_b_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-B"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-B"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-B"].append('MM_'+'_'.join(r_type))
                       
               # HLA-C: 8,9, 30, 31
               r_type = [line[headers.index('r_c_typ1')], line[headers.index('r_c_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_c_typ1')], line[headers.index('d_c_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-C"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-C"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-C"].append('MM_'+'_'.join(r_type))
                       
               # HLA-DRB1: 10,11, 32, 33
               r_type = [line[headers.index('r_drb1_typ1')], line[headers.index('r_drb1_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_drb1_typ1')], line[headers.index('d_drb1_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-DRB1"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-DRB1"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-DRB1"].append('MM_'+'_'.join(r_type))
                       
               # HLA-DQB1: 20,21, 42, 43
               r_type = [line[headers.index('r_dqb1_typ1')], line[headers.index('r_dqb1_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_dqb1_typ1')], line[headers.index('d_dqb1_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-DQB1"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-DQB1"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-DQB1"].append('MM_'+'_'.join(r_type))
                       
               # HLA-DPB1: 24,25, 46, 47
               r_type = [line[headers.index('r_dpb1_typ1')], line[headers.index('r_dpb1_typ2')]]
               r_type.sort()
               d_type = [line[headers.index('d_dpb1_typ1')], line[headers.index('d_dpb1_typ2')]]
               d_type.sort()
               if r_type == d_type:
                   caseIDs["HLA-DPB1"].append('_'.join(r_type))
               else:
                   if (twoFieldHLA(r_type) == twoFieldHLA(d_type)):
                       r_type.extend(d_type)
                       caseIDs["HLA-DPB1"].append('_'.join(r_type))
                   else: 
                       r_type.extend(d_type)
                       caseIDs["HLA-DPB1"].append('MM_'+'_'.join(r_type))

    return(caseIDs)
    
    
    
    
    
    
    
    