#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:00:30 2018

@author: hhuang2
"""

import requests, sys

# http://rest.ensembl.org/#Variation
# serverGR38 = "https://rest.ensembl.org"
rsID = 'rs12157537'
serverGR37 = "http://rest.ensembl.org"
ext = "/variation/homo_sapiens/" + rsID + "?"
 
r = requests.get(serverGR37 + ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))


serverGR38 = "http://rest.ensembl.org"
ext_Clinvar = "/vep/human/region/2:197402635:197402635/T?filter=clinvar"
 
r = requests.get(serverGR38 + ext_Clinvar, headers={ "Content-Type" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))


###
serverGR38 = "https://rest.ensembl.org"
ext_region = "/Homo_sapiens/Gene/Variation_Gene/Table?r=13:32889611-32973805"

r = requests.get(serverGR38 + ext_region, headers={ "Content-Type" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
decoded = r.content
print(repr(decoded))

####
server = "http://rest.ensembl.org"
ext = "/info/variation/consequence_types?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))

######
server = "https://rest.ensembl.org"
ext = "/info/biotypes/homo_sapiens?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))

####### sequence
server = "https://rest.ensembl.org"
ext = "/sequence/region/human/X:1000000..1000100:1?expand_5prime=60;expand_3prime=60"
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)

#
server = "https://rest.ensembl.org"
ext = "/sequence/region/human/X:1000000..1000100:1?mask=soft"
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)

#
server = "https://rest.ensembl.org"
ext = "/sequence/region/human/ABBA01004489.1:1..100?coord_system=seqlevel"
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)

##########
server = "http://rest.ensembl.org"
ext = "/vep/human/region/9:22125503-22125502:1/A?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))
 