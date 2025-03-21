#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 250116

This script takes a tsv file from interproscan as input, 
and outputs a file with a gene ID per line followed by the
location of annotated Pfam domains in the format needed for 
iToL to add to a phylogeny using the 
"dataset_protein_domains_template.txt" file.

"""

__author__ = 'djennah'

#################################################################################################
##Modules##
import argparse

#################################################################################################
#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='parse_interpro_to_iToL',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='{Parse interproscan output to retrieve pfam domains and format for iToL.}')

parser.add_argument('-interpro', '--interpro', required=True, help = 'Path \
	to the input interproscan tsv output file.')

parser.add_argument('-itol', '--itol', required=True, help = 'Gene IDs and \
        pfam information in the format needed for iToL.')

args = parser.parse_args()

#################################################################################################
##Functions##

def get_pfams(interpro):
    with open(interpro, 'r') as f:
        pfam_list = []
        for line in f:
            if line.split("\t")[3] == "Pfam":
                pfam_list.append(line.split("\t")[4])
        return pfam_list

def get_length(interpro):
    with open(interpro, 'r') as f:
        length_dict = {}
        for line in f:
            if line.split("\t")[3] == "Pfam":
                prot_ID=line.split("\t")[0]
                length=line.split("\t")[2]
                length_dict[prot_ID] = length
        return length_dict
    
def parse_interpro(interpro, pfam_dict):
    with open(interpro, 'r') as f:
        interpro_dict = {}
        for line in f:
            if line.split("\t")[3] == "Pfam":
                prot_ID=line.split("\t")[0]
                pfam=line.split("\t")[4]
                pfam_desc=line.split("\t")[4] + " - " + line.split("\t")[5]
                start=line.split("\t")[6]
                stop=line.split("\t")[7]
                shape = pfam_dict[pfam][0]
                colour = pfam_dict[pfam][1]
                entry = '|'.join([shape, start, stop, colour, pfam_desc])
                if prot_ID in interpro_dict:
                    interpro_dict[prot_ID].append(entry) 
                else:
                     interpro_dict[prot_ID] = [entry]
        return interpro_dict

def output_itol(interpro_dict, length_dict, itol):
    for prot_ID in interpro_dict.keys():
        domains = ','.join(interpro_dict[prot_ID])
        itol.write("%s,%s,%s\n" % (prot_ID[:243], length_dict[prot_ID], domains))

#SHAPE|START|END|COLOR|LABEL

#9606,1200,RE|100|150|#ff0000|SH2,EL|400|500|#0000ff|SH3,OC|700|900|#00ff00|PH
    
#################################################################################################
##Implementation##

interpro = args.interpro

#Get pfams and associate with shapes and colours
pfam_list = get_pfams(interpro)
pfam_set = (set(pfam_list))

shapes = ["RE", "HH", "HV", "EL", "DI", "TR", "TL", "PL", "PR", "PU", "PD", "OC", "GP"]
colours = ["#5b859e", "#1e395f", "#75884b", "#1e5a46", "#df8d71", "#af4f2f", "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"]

pfam_dict = {}
num=0
for pfam in pfam_set:
    pfam_dict[pfam] = [shapes[num], colours[num]]
    num+=1
print(pfam_dict)

length_dict = get_length(interpro)
print(length_dict)

#Get pfam domains
interpro_dict = parse_interpro(interpro, pfam_dict)
print(interpro_dict)

#Output file
itol = open(args.itol,'w')
output_itol(interpro_dict, length_dict, itol)