#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 20:38:14 2018

@author: xies
"""

import numpy as np
from Bio import SeqIO, AlignIO, motifs, Alphabet
from Bio.Seq import Seq
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import os
import csv

# ----

# Load Hsap proteome
# Can't use motifs pssm calculate function but can use it as a dictionary
filename = '/data/cycd_targets/hsap_proteome/hsap_proteome.fasta'
ts = SeqIO.parse(filename,'fasta')
hsap_prot  = {}
scores = {}
for rec in ts:
    print rec.name
    fixedID = clean_uniprot_entry(rec.id)
    hsap_prot[fixedID] = rec

# Load Hsap proteome scores
# @todo figure this out later... using either pickle or something like sp savemat
hsap_scores = bg_scores

# Rank targets
max_scores = {rec: np.max(score) for rec, score in hsap_scores.items()}
sorted_scores = sorted(max_scores.iteritems(), key=lambda (k,v): (v,k), reverse=True)
sorted_targets = []
for (k,v) in sorted_scores:
    sorted_targets.append(hsap_prot[k])

# Select one ID to visualize:
idOI = 'P18850'
plt.plot(hsap_scores[idOI])
plt.title(hsap_prot[idOI].name)

# --- Filter for hits above threshold ---
# Additionally filter for unbroken helicity within region
# Filter for presence of at least 1 [ST]P size in protein
threshold = 20
hit_seqs = []
for rec in hsap_scores:
    pos = np.where(hsap_scores[rec] > threshold)[0]
    if pos.size > 0:
        # Check for at least one [S/T]P site (CDK phosphosite)
        if re.search('[ST]P',hsap_prot[rec].seq.tostring()) != None:
            for p in pos:
                # Passed filter; record the seq, pssm scores, and gene name/desc
                s = hsap_prot[rec].seq[p : pssm.length + p]
                score = hsap_scores[rec][p]
                name = hsap_prot[rec].name
                desc = hsap_prot[rec].description
                hit_seqs.append((rec,p,s.tostring(),score,name,desc))

# --- Write scores to a .csv via dataframes
df = pd.DataFrame(hit_seqs,columns=['Entry','Position','Sequence','PSSM score','Name','Description'])
df = df.sort_values('PSSM score',ascending=False)
df.to_csv('/data/cycd_targets/hsap_proteome/hsap_hits>20.csv',sep='\t')

# --- Stringently filter only for F/L-L-R
flr = df[df['Sequence'].apply(flr_filter)]
flr.to_csv( os.path.join(os.path.dirname(filename),'flr.csv'),sep='\t')

def flr_filter(seq):
    s = np.array( list(seq))
    return (s[0] == 'F' or s[0] == 'L') and s[4] == 'L' and s[11] == 'R'
    
