#!/bin/bash
clade=$1
#[ -s GTDB/$clade/profile.txt ] || 
scripts/summarize-phylogenetic-profile.py -id GTDB/$clade/single-copy-ortholog-gene-contexts-pfam -cn GTDB/$clade/orthologs/Results_Aug12/Orthogroups/Orthogroups.GeneCount.tsv --profile GTDB/$clade/profile.txt
for l in upstream downstream both;do
  for d in distinct identical both;do
    #[ -s GTDB/$clade/trait.${l}.${d}.nex ] || 
    scripts/prepare-trait-data.py -i GTDB/$clade/profile.txt -n GTDB/$clade/name.${l}.${d}.txt -o  GTDB/$clade/trait.${l}.${d}.nex -l $l -d $d >  GTDB/$clade/name.${l}.${d}.log 2>&1
    echo GTDB/$clade/trait.${l}.${d}.nex
done
done
