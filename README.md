# tree-utils
Phylogeny utility scripts

## Sequence files, alignments, and trimming

### "remove_seqs_alignment_length.py"

For use in removing short sequences from a multiple sequence alignment (trimmed or untrimmed) that are shorter than a given percentage of the total alignment length.

```
python remove_seqs_alignment_length.py \
  -align "input_alignment" \
  -perc "percentage" \
  -out "output_alignment"
```

### "add_ncbi_taxonomy.py"

For use in adding NCBI taxonomy to fasta sequence headers including NCBI accessions.

First step is to retrieve taxonomy as follows using blastdbcmd from the NCBI BLAST+ package (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and the latest "fullnamelineage.dmp" file from NCBI Taxonomy new_taxdump (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/):

```
while read acc ; do
  taxid=$(blastdbcmd -db nr -dbtype prot -entry $acc -target_only -outfmt '%T') ;
  tax=$(awk -v taxid="$taxid" -F '\t' '$1==taxid {print $5 $3}' fullnamelineage.dmp | sed 's/; /-/g' | sed 's/ /_/g') ;
  echo "$acc $taxid $tax" | tr " " "\t" >> "taxonomy.tsv"
done <"accessions.list"
```

- Compile an "accessions.list" with all NCBI accessions that you want to annotate in your fasta file
- Change "prot" if using nucleotide sequences

Then run the script as follows:

```
python add_ncbi_taxonomy.py \
  -prot "taxonomy.tsv" \
  -faa "sequences.faa" \
  -out "sequences_new-header.faa"
```

## Phylogeny visualization

### "make_marker_tree_pdfs.py"

For use in formatting phylogenies, colouring taxa according to a mapping file, and outputting them as pdfs.

```
  python make_marker_tree_pdfs.py \
    -tree "treefile.tree" \
    -annotation "tree name/annotation" \
    -map "taxa_map.tsv" \
    -outgroup "outgroup.list" \
    -output "treefile.pdf"
```

- For the "taxa_map.tsv", the first column corresponds to the group, the second column to the sequence header, and the third column to the colour
- For the "outgroup.list", list the sequence header for all outgroup sequences that will be used for rooting

### "parse_interpro_to_iToL.py"

For use in parsing InterProScan (https://github.com/ebi-pf-team/interproscan) output files to extract Pfam domains for annotating trees in iToL (https://itol.embl.de/). Append the output to the dataset_protein_domains_template.txt file from iToL (https://itol.embl.de/help.cgi#domains) and drag onto phylogeny.

```
python parse_interpro_to_iToL.py \
  -interpro "interproscan.tsv"  \
  -itol "Pfams_itol_format.txt"
```

- Colours and shapes in script are limited, may need to adjust colours and shapes if needed
