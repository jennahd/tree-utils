#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 220130

This script removes short sequences from a multiple sequence alignment.
"""

__author__ = 'djennah'

#################################################################################################
##Modules##
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#################################################################################################
#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='remove_seqs_alignment_length.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Remove short seqs from an alignment.')

parser.add_argument('-align', '--align', required=True, help = 'Path \
	to the input alignment file.')

parser.add_argument('-perc', '--perc', required=True, help = 'Percentage of alignment cutoff.')

parser.add_argument('-out', '--out', required=True, help = 'Path \
        to the output alignment file with short seqs removed.')

args = parser.parse_args()

#################################################################################################
##Functions##
def read_alignment(alignment):
    alignment = AlignIO.read(open(alignment), "fasta")
    alignment_length = alignment.get_alignment_length()
    print("Alignment length: %i" % alignment_length)
    return alignment, alignment_length

def find_seqs_that_pass_length_cutoff(alignment, perc, alignment_length):
    passed_seq_identifiers = []
    number_filled_needed = int(round(float(alignment_length)*float(perc)/100))
    print("Number of filled positions (not gaps) needed to pass length cutoff: %s" \
    % number_filled_needed)
    print(len(alignment))
    for seq in alignment:
        ID = seq.id.split("|")[0:1][0]
        position = 0
        number_gaps = 0
        number_filled = 0
        while position < alignment_length:
            if seq[position] == "-":
                number_gaps += 1
            else:
                number_filled +=1
            position += 1
        if number_filled >= number_filled_needed:
            passed_seq_identifiers.append(seq.id)
    print(" %s sequences passed the gap length filter" % len(passed_seq_identifiers))
    return passed_seq_identifiers

def output_list_passed_seqs(alignment, passed_seq_identifiers, out):
    out = open(out, "w")
    seq_list = []
    for seq in alignment:
        if seq.id in passed_seq_identifiers:
            seq_list.append(seq)
    SeqIO.write(seq_list, out, "fasta")
    out.close()

########################################################################################################
##Implementation##
alignment = args.align
out = args.out
perc = args.perc

#parse alginment and get alignment length
alignment, alignment_length = read_alignment(alignment)

#find sequences that pass the length cutoff
passed_seq_identifiers = find_seqs_that_pass_length_cutoff(alignment, perc, alignment_length)

#output passed sequeces in a fasta file
output_list_passed_seqs(alignment, passed_seq_identifiers,out)
