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

###### 
# random HML file for 1k genome
import xml.etree.ElementTree as ET
from utils import coreFunctions as cf
from xml.dom import minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def save_hml(prettify_elem, fp):
    
    with open(fp, "w") as f:
        f.write(prettify_elem)
    

#ClassI_template_fp = '../../2018_CDSW_test/ClassI_histo_template.hml'
#ClassII_template_fp = '../../2018_CDSW_test/ClassII_histo_template.hml'

#tree = ET.parse(ClassII_template_fp)
#root = tree.getroot()
#for child in root:
#    print(child.tag, child.attrib)
#    for grandchild in child:
#        print(grandchild.tag, grandchild.attrib)


# read the previous table
HLA_table = pd.read_table(data_fp+'Mockup_igsr_phase1_PairedIDs_HLAtyping-outcome.tsv')
ClassI = ['A', 'B', 'C']
ClassII = ['DRB1', 'DQB1']


### building HML
# 
ROOT = ET.Element('hml',
                  {'project-name': 'GDR_HML_test',
                   'version': '1.0.1',
                   '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation': 'http://schemas.nmdp.org/spec/hml/1.0.1 http://schemas.nmdp.org/spec/hml/1.0.1/hml-1.0.1.xsd'})
ROOT_p1 = ET.SubElement(ROOT,
                        'property', 
                        {'name': 'GeneratedBy', 
                         'value': 'org.nmdp.b2b.hml.service.HMLGatewayService'})

ROOT_p2 = ET.SubElement(ROOT,
                        'property', 
                        {'name': 'MessageReceived', 
                         'value': '2018-04-01 00:00:12.345'})

ROOT_p3 = ET.SubElement(ROOT,
                        'property', 
                        {'name': 'SampleCount', 
                         'value': str(len(HLA_table)*2)}) 

ROOT_p4 = ET.SubElement(ROOT,
                        'property', 
                        {'name': 'UUID', 
                         'value': '123abc-d4e5-6f7g-8h9i-9nmd99u9ucibmtr'})
    
ROOT_p5 = ET.SubElement(ROOT,
                        'hmlid', 
                        {'extension': '999999-2018-4-1', 
                         'root': '999'})
    
ROOT_p6 = ET.SubElement(ROOT,
                        'reporting-center', 
                        {'reporting-center-id': '999'})

## class I
np.random.seed(1)
randINDEX = np.random.permutation(range(len(HLA_table)*2))
counter = 0
Intron2 = 'GTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCTCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACCCCGAGACCCTTGCCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCCGGACGGGGCTCGGGGGACTGGGCTGACCGTGGGGTCGGGGCCAG'
for ind in randINDEX:
    #counter = counter + 1
    if ind < len(HLA_table): # donor
        sampleID = HLA_table['DID'][ind]
        ID_type = 'd'
    else: # recipient
        ind = ind - len(HLA_table)
        sampleID = HLA_table['RID'][ind]
        ID_type = 'r'
    
    Sample = ET.SubElement(ROOT,
                           'sample', 
                           {'center-code': '888', 
                            'id': sampleID})
    
    Sample_cm = ET.SubElement(Sample, 'collection-method')
    Sample_cm.text = 'Whole Blood'
    
    for locus in ClassI:
        counter = counter + 1
        Sample_typing = ET.SubElement(Sample, 'typing',
                                      {'date': '2018-04-01', 
                                       'gene-family': 'HLA'})
        
        Sample_typing_allele = ET.SubElement(Sample_typing, 'allele-assignment',
                                             {'allele-db': 'IMGT/HLA', 
                                              'allele-version': '3.31.0', 
                                              'date': '2018-04-01'})
        tp1 = HLA_table[ID_type+'_'+locus+'_type1_gl'][ind]
        tp2 = HLA_table[ID_type+'_'+locus+'_type2_gl'][ind]
        gl = tp1 +'+'+ tp2
        
        Sample_typing_allele_gl = ET.SubElement(Sample_typing_allele, 'glstring')
        Sample_typing_allele_gl.text = gl
        
        Sample_typing_method = ET.SubElement(Sample_typing, 'typing-method')
        Sample_typing_method_sbt_ngs = ET.SubElement(Sample_typing_method, 'sbt-ngs',
                                              {'locus': 'HLA-'+locus,
                                               'test-id': 'PACBIOSequel',
                                               'test-id-source': 'Histogenetics'})
        
        Sample_typing_consensus = ET.SubElement(Sample_typing, 'consensus-sequence',
                                                {'date': '2018-04-01'})
        
        if tp1 == tp2: # homozygous
            Sample_typing_consensus_ref = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
    
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref,
                                                               'reference-sequence',
                                                               {'id': 'ref'+str(counter), 
                                                                'name': tp1})
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'true', 
                                                                'phase-set': '1', 
                                                                'reference-sequence-id': 'ref'+str(counter)})
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            unaligned_seq = cf.readIMGTsql(tp1, field = 'UnalignedGenomSeq')
            if len(unaligned_seq[0]) == 0:
                seqs = cf.readIMGTsql(tp1, field = 'Exon2,Exon3')
                seqs = Intron2.join(seqs)
                seqs = 'ATGCGTATA' + seqs + 'TTTACGG'
            else:
                seqs = unaligned_seq[0].replace('|','')
                seqs = seqs.replace('*', '')
                
            Sample_typing_consensus_ref_consSeq_seq.text = seqs
         
            
        else: # heterozygous
            # PS-1
            ref_ps1 = 'ref'+str(counter)
            Sample_typing_consensus_ref_ps1 = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref_ps1,
                                                               'reference-sequence',
                                                               {'id': ref_ps1, 
                                                                'name': tp1})
    
            counter = counter + 1
            ref_ps2 = 'ref'+str(counter)
            Sample_typing_consensus_ref_ps2 = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref_ps2,
                                                               'reference-sequence',
                                                               {'id': ref_ps2, 
                                                                'name': tp2})
            # PS1
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'true', 
                                                                'phase-set': '1', 
                                                                'reference-sequence-id': ref_ps1})
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            unaligned_seq = cf.readIMGTsql(tp1, field = 'UnalignedGenomSeq')
            if len(unaligned_seq[0]) == 0:
                seqs = cf.readIMGTsql(tp1, field = 'Exon2,Exon3')
                seqs = Intron2.join(seqs)
                seqs = 'ATGCGTATA' + seqs + 'TTTACGG'
            else:
                seqs = unaligned_seq[0].replace('|','')
                seqs = seqs.replace('*', '')
                
            Sample_typing_consensus_ref_consSeq_seq.text = seqs
            
            # PS-2
            
    
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'true', 
                                                                'phase-set': '2', 
                                                                'reference-sequence-id': ref_ps2})
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            unaligned_seq = cf.readIMGTsql(tp2, field = 'UnalignedGenomSeq')
            if len(unaligned_seq[0]) == 0:
                seqs = cf.readIMGTsql(tp2, field = 'Exon2,Exon3')
                seqs = Intron2.join(seqs)
                seqs = 'ATGCGTATA' + seqs + 'TTTACGG'
            else:
                seqs = unaligned_seq[0].replace('|','')
                seqs = seqs.replace('*', '')
                
            Sample_typing_consensus_ref_consSeq_seq.text = seqs
      
                
save_hml(prettify(ROOT), 'ClassI.hml')


#############################################################
### building HML for Class II
# 
ROOT2 = ET.Element('hml',
                  {'project-name': 'GDR_HML_test',
                   'version': '0.0.1',
                   '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation': 'http://schemas.nmdp.org/spec/hml/1.0.1 http://schemas.nmdp.org/spec/hml/1.0.1/hml-1.0.1.xsd'})
ROOT2_p1 = ET.SubElement(ROOT2,
                        'property', 
                        {'name': 'GeneratedBy', 
                         'value': 'org.nmdp.b2b.hml.service.HMLGatewayService'})

ROOT2_p2 = ET.SubElement(ROOT2,
                        'property', 
                        {'name': 'MessageReceived', 
                         'value': '2018-04-01 00:00:12.345'})

ROOT2_p3 = ET.SubElement(ROOT2,
                        'property', 
                        {'name': 'SampleCount', 
                         'value': str(len(HLA_table)*2)}) 

ROOT2_p4 = ET.SubElement(ROOT2,
                        'property', 
                        {'name': 'UUID', 
                         'value': '123abc-d4e5-6f7g-8h9i-9nmd99u9ucibmtr'})
    
ROOT2_p5 = ET.SubElement(ROOT2,
                        'hmlid', 
                        {'extension': '999999-2018-4-1', 
                         'root': '999'})
    
ROOT2_p6 = ET.SubElement(ROOT2,
                        'reporting-center', 
                        {'reporting-center-id': '999'})

## class II
np.random.seed(1)
randINDEX = np.random.permutation(range(len(HLA_table)*2))
counter = 0
for ind in randINDEX:
    #counter = counter + 1
    if ind < len(HLA_table): # donor
        sampleID = HLA_table['DID'][ind]
        ID_type = 'd'
    else: # recipient
        ind = ind - len(HLA_table)
        sampleID = HLA_table['RID'][ind]
        ID_type = 'r'
    
    Sample = ET.SubElement(ROOT2,
                           'sample', 
                           {'center-code': '888', 
                            'id': sampleID})
    
    Sample_cm = ET.SubElement(Sample, 'collection-method')
    Sample_cm.text = 'Whole Blood'
    
    for locus in ClassII:
        counter = counter + 1
        Sample_typing = ET.SubElement(Sample, 'typing',
                                      {'date': '2018-04-01', 
                                       'gene-family': 'HLA'})
        
        Sample_typing_allele = ET.SubElement(Sample_typing, 'allele-assignment',
                                             {'allele-db': 'IMGT/HLA', 
                                              'allele-version': '3.31.0', 
                                              'date': '2018-04-01'})
        tp1 = HLA_table[ID_type+'_'+locus+'_type1_gl'][ind]
        tp2 = HLA_table[ID_type+'_'+locus+'_type2_gl'][ind]
        gl = tp1 +'+'+ tp2
        Sample_typing_allele_hap1 = ET.SubElement(Sample_typing_allele, 'haploid',
                                                 {'locus': 'HLA-'+locus, 
                                                  'method': 'DNA', 
                                                  'type': tp1.replace(locus+'*', '')})
        Sample_typing_allele_hap1 = ET.SubElement(Sample_typing_allele, 'haploid',
                                                 {'locus': 'HLA-'+locus, 
                                                  'method': 'DNA', 
                                                  'type': tp2.replace(locus+'*', '')})
        
        
        Sample_typing_allele_gl = ET.SubElement(Sample_typing_allele, 'glstring')
        Sample_typing_allele_gl.text = gl
        
        Sample_typing_method = ET.SubElement(Sample_typing, 'typing-method')
        Sample_typing_method_sbt_ngs = ET.SubElement(Sample_typing_method, 'sbt-ngs',
                                              {'locus': 'HLA-'+locus,
                                               'test-id': 'RSII',
                                               'test-id-source': 'Histogenetics'})
        
        Sample_typing_consensus = ET.SubElement(Sample_typing, 'consensus-sequence',
                                                {'date': '2018-04-01'})
        
        if tp1 == tp2: # homozygous
            Sample_typing_consensus_ref = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref,
                                                               'reference-sequence',
                                                               {'id': 'ref'+str(counter), 
                                                                'name': tp1})
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon2', 
                                                                'phase-set': '1', 
                                                                'reference-sequence-id': 'ref'+str(counter)})
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            ARS_seq = cf.readIMGTsql(tp1, field = 'Exon2,Exon3')
            Sample_typing_consensus_ref_consSeq_seq.text = ARS_seq[0]
            
            if len(ARS_seq[1].replace('*', '')) > 0: # if Exon 3 is not empty
                #PS_count += 1
                Sample_typing_consensus_ref_consSeq_3 = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon3',
                                                                'phase-set': '1',
                                                                'reference-sequence-id': 'ref'+str(counter)})
    
                Sample_typing_consensus_ref_consSeq_seq_3 = ET.SubElement(Sample_typing_consensus_ref_consSeq_3,
                                                                    'sequence')
                Sample_typing_consensus_ref_consSeq_seq_3.text = ARS_seq[1].replace('*', '')
         
            
        else: # heterozygous
            # PS-1
            ref_ps1 = 'ref'+str(counter)
            Sample_typing_consensus_ref_ps1 = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref_ps1,
                                                               'reference-sequence',
                                                               {'id': ref_ps1, 
                                                                'name': tp1})
    
            counter = counter + 1
            ref_ps2 = 'ref'+str(counter)
            PS_count = 1
            Sample_typing_consensus_ref_ps2 = ET.SubElement(Sample_typing_consensus, 'reference-database',
                                                    {'availability': 'public', 
                                                     'name': 'IMGT/HLA', 
                                                     'version': '3.31.0'})
            Sample_typing_consensus_ref_refseq = ET.SubElement(Sample_typing_consensus_ref_ps2,
                                                               'reference-sequence',
                                                               {'id': ref_ps2, 
                                                                'name': tp2})
            # PS1
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon2',
                                                                'phase-set': str(PS_count),
                                                                'reference-sequence-id': ref_ps1})
    
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            ARS_seq = cf.readIMGTsql(tp1, field = 'Exon2,Exon3')
            Sample_typing_consensus_ref_consSeq_seq.text = ARS_seq[0]
            
            if len(ARS_seq[1].replace('*', '')) > 0: # if Exon 3 is not empty
                #PS_count += 1
                Sample_typing_consensus_ref_consSeq_3 = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon3',
                                                                'phase-set': str(PS_count),
                                                                'reference-sequence-id': ref_ps1})
    
                Sample_typing_consensus_ref_consSeq_seq_3 = ET.SubElement(Sample_typing_consensus_ref_consSeq_3,
                                                                    'sequence')
                Sample_typing_consensus_ref_consSeq_seq_3.text = ARS_seq[1].replace('*', '')


            # PS-2
            PS_count += 1
            Sample_typing_consensus_ref_consSeq = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon2',
                                                                'phase-set': str(PS_count),
                                                                'reference-sequence-id': ref_ps2})
            Sample_typing_consensus_ref_consSeq_seq = ET.SubElement(Sample_typing_consensus_ref_consSeq,
                                                                    'sequence')
            
            ARS_seq = cf.readIMGTsql(tp2, field = 'Exon2,Exon3')
            Sample_typing_consensus_ref_consSeq_seq.text = ARS_seq[0]
            
            if len(ARS_seq[1].replace('*', '')) > 0: # if Exon 3 is not empty
                #PS_count += 1
                Sample_typing_consensus_ref_consSeq_3 = ET.SubElement(Sample_typing_consensus,
                                                               'consensus-sequence-block',
                                                               {'continuity': 'false',
                                                                'description': 'HLA-'+locus+' Exon3',
                                                                'phase-set': str(PS_count),
                                                                'reference-sequence-id': ref_ps2})
    
                Sample_typing_consensus_ref_consSeq_seq_3 = ET.SubElement(Sample_typing_consensus_ref_consSeq_3,
                                                                    'sequence')
                Sample_typing_consensus_ref_consSeq_seq_3.text = ARS_seq[1].replace('*', '')
      
                
save_hml(prettify(ROOT2), 'ClassII.hml')





