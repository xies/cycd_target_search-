#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 18:23:34 2018

@author: xies
"""

import numpy as np
import scipy as sp
from Bio import SeqIO, AlignIO, motifs, Alphabet
from Bio.Seq import Seq
import pandas as pd
import os

# --- Load target FASTA sequences ---
filename = '/data/database/uniprot-reviewed%3Ayes+AND+proteome%3Aup000005640.fasta'
ts = SeqIO.parse(filename,'fasta')
targets  = {}
scores = {}
for rec in ts:
    print rec.name
    scores[rec.id] = calculate_scores(pssm, rec.seq )
    targets[rec.id] = rec

# --- Load SS pred ----
filename = '/data/cycd_targets/cycd_target_uniprot_wider_individuals/psipred.fasta'
ss = SeqIO.parse(filename,'fasta')
secondary = {}
for rec in ss:
    secondary[rec.id] = rec

    
# Further filter by a threshold of helicity and record all as dataframe
filtered_hits = []
for k in hit_seqs:
    hits = hit_seqs[k]
    hits = zip(hits[::3],hits[1::3],hits[2::3])
    for (s,h,score) in hits:
        helicity = np.array( list(h) )
        I = helicity == 'H'        
        if sum(I) > 8: # If greater than 8 were helical, check if there are breaks in helicity
            hbreak = np.diff(np.cumsum(I))
            if all(hbreak[1:-1] != 0):  # breaks other than ends
                # Grab the sequence metadata
                name = targets[k].name
                description = targets[k].description
                filtered_hits.append( (k,s,h,score,name,description) )
df = pd.DataFrame(filtered_hits,columns=['Entry','Sequence','Helicity','PSSM score','Name','Description'])

df.sort_values('PSSM score',ascending=False).to_csv( os.path.join(os.path.dirname(filename),'hits.csv'),sep='\t')
