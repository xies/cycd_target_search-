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

# Load hand-curated CycD putative target list
filename = '/data/cycd_targets/cycd_target_uniprot.txt'
targetIDs = pd.read_csv(filename)
entries = targetIDs['Entry']

# Fetch and write as FASTA
out_name = '/data/cycd_targets/cycd_target_uniprot.fasta'
upData = uniprot.batch_uniprot_metadata(entries, 'cache')
uniprot.write_fasta(out_name, upData, entries)

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




