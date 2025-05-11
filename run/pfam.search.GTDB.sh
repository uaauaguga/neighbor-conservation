#!/bin/bash
#for clade in $(ls GTDB  );do
for clade in Actinomycetia;do
   # [ -d GTDB/$clade/proteins.pfam ] || 
   sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=10 --nodes=1 --wrap "bash scripts/pfam.search.GTDB.sh $clade "
done
