#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 17:38:57 2018

@author: xies
"""

import numpy as np
import pandas as pd
import os, re, copy
from Bio import SeqIO,AlignIO,Seq

PSIPRED_DIR = '/data/cycd_targets/cycd_target_uniprot_individuals'

seqs = []

for filename in os.listdir(PSIPRED_DIR):
    if filename.endswith('.ss2'):
        print 'Working on ', filename
        fastaname,ext = os.path.splitext(filename)
        
        filename = os.path.join(UNGAPPED_DIR,filename)
        
        #Load PSIPRED VFORMAT in a sane way to extract only relevant info
        df = pd.read_csv(filename,header=0, delim_whitespace=True,skiprows=0,
                         names=['position','AA','prediction'],usecols=[0,1,2], index_col=0)
        
        seq = Seq.Seq(''.join(df.prediction))
        seq.id = os.path.splittext()
        # Load the same GAPPED sequence file
#        seq = SeqIO.read(os.path.join(GAPPED_DIR,fastaname),'fasta')
#        seq_ss = copy.deepcopy(seq) # deepcopy to get a hard copy
        
        # Find parts of sequence that's not the gapped character
#        s = seq.seq.tostring()
#        s_array = np.array(list(s))
#        I = s_array != '-'
#        s_array[I] = df.prediction
#        seq_ss.seq = Seq.Seq(''.join(s_array))
#        seq_ss.id = os.path.splitext(fastaname)[0]
#        seqs.append(seq_ss)
#        
    else:
        continue

SeqIO.write(seqs,os.path.join( GAPPED_DIR,'psipred.fasta' ),'fasta')
