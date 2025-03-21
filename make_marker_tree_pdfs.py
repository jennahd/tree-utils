#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 221017

This script performs basic formatting of newick format trees and prints them as
a pdf. If you have any defined labels in a mapping file it will colour those branches,
and it will root by a given outgroup.

"""
#Modules
from ete3 import *
import argparse
import os
from collections import Counter

###Command-line usage with argparser

parser = argparse.ArgumentParser(prog='make_marker_tree_pdfs.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Print newick tree file as pdf with taxa colour labelled.')

parser.add_argument('-tree', '--tree', required=True, help = 'File containing \
	tree in Newick format')

parser.add_argument('-annotation', '--annotation', required=True, help = 'Title for tree \
	to be printed in pdf with gene and corresponding annotation')

parser.add_argument('-map', '--map', required=True, help = 'Colour mapping file \
	with each taxon prefix and corresponding taxon colour')

parser.add_argument('-outgroup', '--outgroup', required=True, help = 'File with \
	list of outgroup taxa.')

parser.add_argument('-output', '--output', required=True, help = 'File for tree \
	to be written to pdf.')

args = parser.parse_args()

print(args)

##Check for files and arguments
if not os.path.exists(args.tree):
	raise IOError("This file does not exist")

##############################################################################
#Functions

def file_to_list(file):
    list = []
    with open(file, 'r') as f:
        for line in f:
            list.append(line.strip("\n"))
    return list

def parse_mapping_file(map):
    colour_dict = {}
    names = []
    for line in map:
        line = line.strip()
        group, name, colour = line.split("\t")[0:3]
        colour_dict[name] = group, colour
        names.append(name)
    return colour_dict, names

def initialize_tree(tree):
    tree = Tree(args.tree, format=0)
    tree.ladderize(direction=1)
    leaves = tree.get_leaf_names()
    return tree, leaves

def find_outgroup_seqs(leaves, outgroup_names):
    outgroup_leaves = []
    for leaf in leaves:
        name = "@".join(leaf.split('@', 2)[:2])
        if name in outgroup_names:
            outgroup_leaves.append(leaf)
    return(outgroup_leaves)

def reroot_tree(tree, outgroup_leaves):
    monophyly = tree.check_monophyly(values=outgroup_leaves, target_attr="name", ignore_missing=True, unrooted=True)
    monophyly = str(monophyly[0])
    if monophyly == 'True':
        print("All outgroup taxa are monophyletic")
        outgroup = tree.get_common_ancestor(outgroup_leaves)
        tree.set_outgroup(outgroup)
    else:
        print("Not all outgroup taxa are monophyletic")
        outgroup_mid = tree.get_midpoint_outgroup()
        tree.set_outgroup(outgroup_mid)
    return tree

def count_dict(tree, names, leaves):
    counts = {}
    for leaf in leaves:
        name = "@".join(leaf.split('@', 2)[:2])
        if name in counts.keys():
            print("duplicate taxon: %s" %(name))
            counts[name][0].append(leaf)
            counts[name][1] += 1
        else:
            counts[name] = [[leaf], 1]
    return counts

def duplicates(counts):
    duplicates = []
    for name in counts.keys():
        if counts[name][1] >= 2:
            duplicates = duplicates + counts[name][0]
    return duplicates

def label_duplicates(duplicates, tree):
    style = NodeStyle()
    style["fgcolor"] = "#000000"
    style["size"] = 0
    style["hz_line_width"] = 2
    style["hz_line_color"] = "#000000"
    style["vt_line_width"] = 2
    style["vt_line_color"] = "#000000"
    dstyle = NodeStyle()
    dstyle["fgcolor"] = "#FF0000"
    dstyle["size"] = 10
    dstyle["hz_line_width"] = 5
    dstyle["hz_line_color"] = "#FF0000"
    for node in tree.traverse():
        if node.is_leaf():
            if node.name in duplicates:
                node.set_style(dstyle)
            else:
                node.set_style(style)
        else:
            node.set_style(style)
    return tree

def add_tree_colours(tree, colour_dict):
    for node in tree.traverse():
        if node.is_leaf():
            name = "@".join(node.name.split('@', 2)[:2])
            name_face = TextFace(node.name, fgcolor=colour_dict[name][1])
            name_face2 = TextFace(colour_dict[name][0] + "@", fgcolor= colour_dict[name][1])
            node.add_face(name_face, column=1, position='branch-right')
            node.add_face(name_face2, column=0, position='branch-right')
    return tree

def output_tree(tree, output, annotation):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.title.add_face(TextFace(annotation,fsize=12,bold=True),0)
    ts.scale = 200
    tree.render(args.output,tree_style=ts)

##############################################################################
#Implementation

#Annotation
annotation = "%s" %(args.annotation)
print(annotation)

#Make lists from files
outgroup_names = file_to_list(args.outgroup)
map = file_to_list(args.map)

#Parse mapping file
colour_dict, names = parse_mapping_file(map)

#Initialize tree
tree, leaves = initialize_tree(args.tree)

#Identify outgroup seqs
outgroup_leaves = find_outgroup_seqs(leaves, outgroup_names)

#Reroot with outgroup
tree = reroot_tree(tree, outgroup_leaves)

#Count duplicate taxa
counts = count_dict(tree, names, leaves)
duplicates = duplicates(counts)

#Colour duplicates in tree
tree = label_duplicates(duplicates, tree)

#Colour taxa in tree
tree = add_tree_colours(tree, colour_dict)

#Output tree pdf
output_tree(tree, args.output, annotation)
