#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 18:12:37 2018

@author: xies
"""

import numpy as np
import scipy as sp
import uniprot
from Bio import SeqIO, AlignIO, motifs, Alphabet
from Bio.Seq import Seq
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import os

# Load excel file of all putative targets
filename = '/data/cycd_targets/1-s2.0-S1535610811003606-mmc2.xls'
df = pd.read_excel(filename)
up = df['UniProKB/SWISS-PROT'].str.split(', ',1).tolist()[1:] # split the Uniprot entry field
up = pd.DataFrame(np.vstack( up[1:] ), columns=['Entry name','Entry']) # Extract into new dataframe
entries = up['Entry'].str.strip(' ') # Make sure UniProt entry ID has no spaces

# Load hand-curated CycD putative target list
filename = '/data/cycd_targets/cycd_target_uniprot.txt'
targetIDs = pd.read_csv(filename)
already_seen = pd.concat((targetIDs['Entry'],entries))

# Load hit list from PSSM
filename = '/data/cycd_targets/hsap_proteome/hsap_hits>20.csv'
targetIDs = pd.read_csv(filename,sep='\t')
entries = targetIDs['Entry']

# Do a merge to see what's not already seen in the hand-curated list
merged = set(entries) - set(already_seen)

# Fetch and write as FASTA
out_name = '/data/cycd_targets/hsap_hits>20.fasta'
upData = uniprot.batch_uniprot_metadata(merged, 'cache')
uniprot.write_fasta(out_name, upData, merged)

split_fastas(out_name)


PSIPRED_DIR = '/data/cycd_targets/cycd_target_uniprot_individuals'
seqs = []

for filename in os.listdir(PSIPRED_DIR):
    if filename.endswith('.ss2'):
        print 'Working on ', filename
        
        #Load PSIPRED VFORMAT in a sane way to extract only relevant info
        df = pd.read_csv(os.path.join(PSIPRED_DIR,filename),header=0, delim_whitespace=True,skiprows=0,
                         names=['position','AA','prediction'],usecols=[0,1,2], index_col=0)
        fastaname,ext = os.path.splitext(filename)
        s = SeqIO.read(os.path.join(PSIPRED_DIR,fastaname),'fasta')
        
        s.seq = Seq(''.join(df.prediction))
        seqs.append(s)
        
    else:
        continue

SeqIO.write(seqs,os.path.join( PSIPRED_DIR,'psipred.fasta' ),'fasta')




