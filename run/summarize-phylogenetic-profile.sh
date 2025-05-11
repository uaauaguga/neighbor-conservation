#!/bin/bash
#for clade in $(ls GTDB);do
for clade in Actinomycetia;do
  sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=1 --nodes=1 --wrap "bash scripts/summarize-phylogenetic-profile.sh $clade" 
done
