#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 18:38:08 2018

@author: xies
"""

import numpy as np
from Bio import SeqIO, AlignIO, motifs, Alphabet
from Bio.Seq import Seq
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import os

# -----

# Load RB helix sequences
filename = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/manually_filtered/mafft_trimmed.fasta'
helix = AlignIO.read(filename,'fasta')
helix = [rec.seq for rec in helix]

# Massage sequence to fit IUPAC protein alphabet
instances = []
for rec in helix:
    s = np.array( list(rec.tostring()) )
    s[ s == '-' ] = 'X'
    s = s[23:35]
    rec = Seq(''.join(s),Alphabet.IUPAC.ExtendedIUPACProtein)
#    rec.alphabet = Alphabet.IUPAC.ExtendedIUPACProtein
    instances.append(rec)

# Calculate counts matrix, normalize and log
m = motifs.create(instances)
pwm = m.counts.normalize(pseudocounts=.5)
pssm = pwm.log_odds()

# Write PSSM to .csv via dataframes
filename = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/manually_filtered/mafft_trimmed_pssm.csv'
df = pd.DataFrame(pssm)
df.to_csv(filename)

# ----

# Load cycD target sequences
# Can't use motifs pssm calculate function but can use it as a dictionary
filename = '/data/cycd_targets/cycd_target_uniprot_wider.fasta'
ts = SeqIO.parse(filename,'fasta')
targets  = {}
scores = {}
for rec in ts:
    print rec.name
    fixedID = clean_uniprot_entry(rec.id)
    scores[fixedID] = calculate_scores(pssm, rec.seq )
    targets[fixedID] = rec

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

# Estimate PSSM background with randomized sequences

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