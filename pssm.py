#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 18:38:08 2018

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
from random import shuffle

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

# ----

# Load cycD target sequences
# Can't use motifs pssm calculate function but can use it as a dictionary
filename = '/data/database/uniprot-reviewed%3Ayes+AND+proteome%3Aup000005640.fasta'
ts = SeqIO.parse(filename,'fasta')
targets  = {}
scores = {}
for rec in ts:
    print rec.name
    scores[rec.id] = calculate_scores(pssm, rec.seq )
    targets[rec.id] = rec


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
# Concatenate all sequences
randseq = []
for rec in targets:
    randseq.append( targets[rec].seq.tostring() )
randseq = ''.join(randseq)
randseq = list(randseq)
shuffle(randseq)
randseq = ''.join(randseq) # in-place
randscores = calculate_scores(pssm,randseq)
mean_randscore = np.mean(randscores)
std_randscore = np.std(randscores)
print 'Mean random score: ', mean_randscore
print 'Std random score: ', std_randscore
print '3.5 std above mean: ', mean_randscore + 3.5 * std_randscore

# --- Load SS pred ----
filename = '/data/cycd_targets/cycd_target_uniprot_wider_individuals/psipred.fasta'
ss = SeqIO.parse(filename,'fasta')
secondary = {}
for rec in ss:
    secondary[rec.id] = rec

# Apply threshold
threshold = 3
hit_seqs = {}
for rec in scores:
    pos = np.where(scores[rec] > threshold)[0]
    if pos.size > 0:
        try:
            hits = []
            for p in pos:
                # Record the seq, the second struct, and pssm scores
                hits.append(targets[rec].seq[p : pssm.length + p].tostring())
                hits.append(secondary[rec].seq[p : pssm.length + p].tostring())
                hits.append(scores[rec][p])
            hit_seqs[rec] = hits
        except KeyError:
            continue
        

# Print to shell (useful for small number of candidates)
for k in hit_seqs:
    print '---'
    print targets[k].description
    print hit_seqs[k]
    
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

# ---
    
def calculate_scores(pssm,s):
    
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
    
    