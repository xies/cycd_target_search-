#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 20:34:44 2018

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

# Load cycD target sequences
# Can't use motifs pssm calculate function but can use it as a dictionary
filename = '/data/cycd_targets/putative_targets.fasta'
ts = SeqIO.parse(filename,'fasta')
targets  = {}
scores = {}
for rec in ts:
    print rec.name
    fixedID = clean_uniprot_entry(rec.id)
    scores[fixedID] = calculate_scores(pssm, rec.seq )
    targets[fixedID] = rec
write_dict_to_csv(scores,'/data/cycd_targets/putative_targets/scores.csv')

# Visualize the distribution of scores
fig,ax = plt.subplots()
for rec in scores:    
    ax.plot(scores[rec])
    
plt.hist( np.hstack(scores.values()) ,1000)

# Rank targets
max_scores = {rec: np.max(score) for rec, score in scores.items()}
sorted_scores = sorted(max_scores.iteritems(), key=lambda (k,v): (v,k), reverse=True)
sorted_targets = []
for (k,v) in sorted_scores:
    sorted_targets.append(targets[k])

# Select one ID to visualize:
idOI = 'P18850'
plt.plot(scores[idOI])
plt.title(targets[idOI].name)

#---- BG estimation
# Estimate PSSM background with randomized sequences
#bg_scores = estimate_random_background(targets,pssm)

# Estimate PSSM background with Hsap genome
filename = '/data/cycd_targets/hsap_proteome/hsap_proteome.fasta'
ts = SeqIO.parse(filename,'fasta')
bg  = {}
bg_scores = {}
for rec in ts:
    print rec.name
    fixedID = clean_uniprot_entry(rec.id)
    bg_scores[fixedID] = calculate_scores(pssm, rec.seq )
    bg[fixedID] = rec
    
# Save background scores for future use
write_dict_to_csv(bg_scores,'/data/cycd_targets/pssm_hsap_background.csv')

bgAll = np.hstack(bg_scores.values())
plt.hist( bgAll ,1000)

print 'Mean random score: ', np.mean(bgAll)
print 'Std random score: ', np.std(bgAll)
print '3.5 std above mean: ', np.mean(bgAll) + 3.5 * np.std(bgAll)

# --- Filter for hits above threshold
threshold = 20
hit_seqs = []
for rec in scores:
    pos = np.where(scores[rec] > threshold)[0]
    if pos.size > 0:
        print rec
        hits = []
        for p in pos:
            # Record the seq, pssm scores, and gene name/desc
            hits.append(rec)
            hits.append(targets[rec].seq[p : pssm.length + p].tostring())
            hits.append(scores[rec][p])
            hits.append(targets[rec].name)
            hits.append(targets[rec].description)
            hit_seqs.append(hits)

# --- Write scores to a .csv via dataframes
df = pd.DataFrame(hit_seqs,columns=['Entry','Sequence','PSSM score','Name','Description'])
df.sort_values('PSSM score',ascending=False).to_csv('/data/cycd_targets/hsap_proteome/hsap_hits>20.csv',sep='\t')

# ---
    
def calculate_scores(pssm,s):
    # Apply PSSM to a specific sequence s
    s = s.upper()
    scores = []

    seqLen = len(s)
    print seqLen
    motifLen = pssm.length
    
    for i in range(seqLen - motifLen + 1):
        score = 0.0
        for position in range(motifLen):
            letter = s[i + position]
            try:
                score += pssm[letter][position]
            except KeyError:
                score = float('nan')
                break
        scores.append(score)
    
    return np.array(scores)
    
def estimate_random_background(targets,pssm):
    # Randomize input sequences and calculate background distribution of scores
    # for a given PSSM
    

    from random import shuffle
    # Concatenate all sequences
    randseq = []
    for rec in targets:
        randseq.append( targets[rec].seq.tostring() )
    randseq = ''.join(randseq)
    randseq = list(randseq)
    # shuffle-in place
    shuffle(randseq)
    randseq = ''.join(randseq) # join back up again
    randscores = calculate_scores(pssm,randseq)
    mean_randscore = np.mean(randscores)
    std_randscore = np.std(randscores)
    print 'Mean random score: ', mean_randscore
    print 'Std random score: ', std_randscore
    print '3.5 std above mean: ', mean_randscore + 3.5 * std_randscore
    
    return randscores
