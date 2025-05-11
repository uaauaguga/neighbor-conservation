#!/bin/bash
clade=$1
mkdir -p GTDB/$clade/single-copy-ortholog-gene-contexts
mkdir -p GTDB/$clade/single-copy-ortholog-gene-contexts-pfam
for bed in $(ls GTDB/$clade/genes );do
  asm_id=${bed%.*}
  echo $asm_id
  #[ -s GTDB/$clade/single-copy-ortholog-gene-contexts/${asm_id}.txt ] || 
   scripts/gather-contexts-GTDB.py -i GTDB/$clade/genes/${asm_id}.bed --output GTDB/$clade/single-copy-ortholog-gene-contexts/${asm_id}.txt --og-annotation GTDB/$clade/orthologs/Results_Aug12/Orthogroups/Orthogroups.tsv --pfam-annotation GTDB/$clade/proteins.pfam/${asm_id}.txt
  #[ -s GTDB/$clade/single-copy-ortholog-gene-contexts-pfam/${asm_id}.txt ] || 
  scripts/summarize-context-pfam.py -i GTDB/$clade/single-copy-ortholog-gene-contexts/${asm_id}.txt -o GTDB/$clade/single-copy-ortholog-gene-contexts-pfam/${asm_id}.txt 
done
echo 
  #[ -s GTDB/$clade/profile.${mif}.txt ] || 
  #scripts/summarize-phylogenetic-profile.py --input-directory GTDB/$clade/single-copy-ortholog-gene-contexts-pfam  --copy-number GTDB/$clade/orthologs/Results_Aug12/Orthogroups/Orthogroups.GeneCount.tsv --profile GTDB/$clade/profile.txt 
