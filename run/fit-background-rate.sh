#!/bin/bash
for clade in $(ls GTDB);do
for l in both;do #upstream downstream both;do
  for d in both;do #distinct identical both;do
    echo GTDB/$clade/trait.${l}.${d}.nex
    sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=4 --nodes=1 --wrap "rb scripts/estimate-rate.free.binary.uniform.Rev --args GTDB/$clade/trait.${l}.${d}.nex GTDB/$clade/bac120.nwk GTDB/$clade/rate.${l}.${d}.uniform.20.txt"
  done
done
done
