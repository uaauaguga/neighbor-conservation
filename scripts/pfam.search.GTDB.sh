#!/bin/bash
clade=$1
mkdir -p GTDB/$clade/proteins.pfam
for faa in $(ls GTDB/$clade/proteins );do
  asm_id=${faa%.*}
  scripts/annotate-proteins-to-pfam.py -q GTDB/$clade/proteins/$faa -o GTDB/$clade/proteins.pfam/${asm_id}.txt -td tmp/$clade > GTDB/$clade/proteins.pfam/${asm_id}.log 2>&1
done
