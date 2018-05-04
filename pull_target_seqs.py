#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 18:12:37 2018

@author: xies
"""

import numpy as np
import scipy as sp
import uniprot
from Bio import SeqIO, AlignIO
import pandas as pd
import seaborn as sb

# Load hand-curated CycD putative target list
filename = '/data/cycd_targets/cycd_target_uniprot.txt'
targetIDs = pd.read_csv(filename)
entries = targetIDs['Entry']

# Fetch and write as FASTA
out_name = '/data/cycd_targets/cycd_target_uniprot.fasta'
upData = uniprot.batch_uniprot_metadata(entries, 'cache')
uniprot.write_fasta(out_name, upData, entries)

split_fastas(out_name)


