#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 18:43:32 2018

@author: xies
"""



def split_fastas(filename):
    #Split a .fasta file with multiple entries into multiple
    #individual .fasta files
    import os
    from Bio import SeqIO
     
    # Get file and directory names for unified output
    indir = os.path.dirname(filename)
    base_filename = os.path.basename(filename)
    base_noext = os.path.splitext(base_filename)[0]
    base_noext = ''.join((base_noext,'_individuals'))
    out_dirname = os.path.join(indir,base_noext)
    # make output directory if it DNE
    try:
        os.stat(out_dirname)
    except:
        os.mkdir(out_dirname)

    seqs = SeqIO.parse(filename,'fasta')
    
    for rec in seqs:
       rec_basename = ''.join((rec.name,'.fasta'))
       rec_filename = os.path.join(out_dirname,rec_basename)
       rec.seq = rec.seq.ungap('-')
       if rec.seq == '':
           continue
       print rec_filename
       with open(rec_filename,'w') as output_handle:
           SeqIO.write(rec,output_handle,'fasta')
           
def clean_uniprot_entry(badentry):
    # Grab the UniProt Entry field from the sp|ENTRY|NAME format of downloaded FASTAs
    s = badentry.split('|')
    if len(s) > 2:
        goodentry = s[1]
    else:
        goodentry = badentry
    return goodentry


def write_dict_to_csv(dict,filename):
    # Save dictionary for future use
    import csv
    with open(filename,'wb') as csv_file:
        writer = csv.writer(csv_file)
        for k,v in dict.items():
            writer.writerow([k,v])

