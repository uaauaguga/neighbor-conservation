#!/bin/bash
for gff in $(ls genome/annotations);do
  asm_id=$(echo $gff | cut -f 1,2 -d '_')
  #scripts/gff2bed.py --gff genome/annotations/$gff --bed genome/bed/${asm_id}.bed --feature CDS --name protein_id 
  sort -k1,1 -k2,2n -o genome/bed/${asm_id}.bed genome/bed/${asm_id}.bed
done
