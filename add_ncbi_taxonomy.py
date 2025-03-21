#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 220120

This script takes input ncbi accessions and taxids and appends taxonomy to fasta sequence files.
"""

__author__ = 'djennah'

#################################################################################################
##Modules##
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

#################################################################################################
#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='add_ncbi_taxonomy.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Adds taxonomy to fasta files based on taxids.')

parser.add_argument('-prot', '--prot', required=True, help = 'Path \
	to the input tsv file with protein accessions and taxonomy.')

parser.add_argument('-faa', '--faa', required=True, help = 'Path \
        to the input faa file to add taxonomy to.')

parser.add_argument('-out', '--out', required=True, help = 'Path \
        to the output faa file.')

args = parser.parse_args()

#################################################################################################
##Functions##
def get_tax(prot):
    with open(prot, 'r') as f:
        tax_dict = {}
        for line in f:
            accession = line.split("\t")[0:1][0]
            try:
                tax = line.split("\t")[2:3][0].strip("\n")
            except:
                tax = "na"
            tax_dict[accession] = tax
        return tax_dict

def replace_fasta_headers(faa, tax_dict, out):
    output_file = open(out, "w")
    taxonomy = "na"
    for seq_record in SeqIO.parse(faa, "fasta"):
        accession = seq_record.description.split(" ",1)[0].replace("NR@","")
        try:
            description = seq_record.description.split(" ",1)[1].split("]",1)[0].replace(" ","_").replace("(","").replace(")","").replace(":","_").replace("\\","_").replace(";","_").replace("/","_").replace("=","_") + "]"
        except:
            try:
                description = seq_record.description.split(" ",1)[1].replace(" ","_").replace("(","").replace(")","").replace(":","_").replace("\\","_").replace(";","_").replace("/","_").replace("=","_")
            except:
                description = "na"
        if accession in tax_dict.keys():
            taxonomy = tax_dict[accession]
        else:
            taxonomy = "na"
        new_seq_header = seq_record.description.split(" ",1)[0] + "@" + description + "@" + taxonomy
        new_seq_header.replace(":","_").replace(";","_").replace("(","_").replace(")","_")
        new_seq_header = new_seq_header
        output_file.write(">%s\n%s\n" % (
        new_seq_header,
        seq_record.seq))

#################################################################################################
##Implementation##
prot = args.prot
faa = args.faa
out = args.out

#get taxids
tax_dict = get_tax(prot)

#Replace fasta headers
replace_fasta_headers(faa, tax_dict, out)
