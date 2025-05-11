#!/bin/bash
#for clade in $(ls GTDB);do
#for clade in Actinomycetia;do
for clade in 8-genomes;do 
  sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=1 --nodes=1 --wrap "bash scripts/gather-contexts-GTDB.sh $clade" 
done
