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
#scores = load...

# --- Load SS pred ----
filename = '/data/cycd_targets/cycd_target_uniprot_wider_individuals/psipred.fasta'
ss = SeqIO.parse(filename,'fasta')
secondary = {}
for rec in ss:
    secondary[rec.id] = rec


# --- Filter for hits above threshold
threshold = 5
hit_seqs = []
for rec in scores:
    pos = np.where(scores[rec] > threshold)[0]
    if pos.size > 0:
        print rec
        for p in pos:
            # Pull the secondary structure prediction
            helicity = np.array(list( secondary[rec][p:pssm.length + p].seq.tostring() ))
            I = helicity == 'H'
            if sum(I) > 9: # Check if greater than 8 are helical
                # Check if there's internal break in helicity
                hbreak = np.diff(np.cumsum(I))
                if all(hbreak[1:-1] != 0):
                    # Passed filter; record the seq, pssm scores, and gene name/desc
                    s = targets[rec].seq[p : pssm.length + p]
                    h = ''.join(helicity)
                    score = scores[rec][p]
                    name = targets[rec].name
                    desc = targets[rec].description
                    hit_seqs.append((rec,s.tostring(),h,score,name,desc))

df = pd.DataFrame(hit_seqs,columns=['Entry','Sequence','Helicity','PSSM score','Name','Description'])

df.sort_values('PSSM score',ascending=False).to_csv( os.path.join(os.path.dirname(filename),'hits.csv'),sep='\t')


