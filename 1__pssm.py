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
import csv

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
