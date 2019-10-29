#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Script to create/recreate peptide file
'''

import pandas as pd
from pyteomics import fasta, parser

#parameters
fastaPath = './sprot_human.fasta' #path to FASTA file
minlength = 6 #minimal peptide length
maxlength = 20 #maximal peptide length
maxMC = 1 #maximum number of missed cleavages

peptides = []
with fasta.read(fastaPath) as f:
    for n,s in f:
        peptides += parser.cleave(s,rule=parser.expasy_rules['trypsin'], \
                                  min_length=minlength, missed_cleavages=maxMC)

peptides = pd.Series([x for x in set(peptides) if len(x) < maxlength], name='sequence')
badAA = peptides['sequence'].str.findall("[^ACDEFGHIKLMNPQRSTVWY]").apply(len) > 0 #remove non-standard AA
peptides.drop(peptides.index[badAA], axis='index', inplace=True)
peptides.drop_duplicates().to_csv('./helps/peptides.csv', index=False) #remove duplicate sequences
print("Number of peptides in the list: {}".format(len(peptides)))
