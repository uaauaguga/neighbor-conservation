#!/bin/bash
for faa in $(ls genome/proteins);do
  asm_id=$(echo $faa | cut -f 1,2 -d '_')
  scripts/annotate-proteins-to-pfam.py -q genome/proteins/$faa -o genome/proteins.pfam/${asm_id}.txt -td tmp > genome/proteins.pfam/${asm_id}.log 2>&1
done
