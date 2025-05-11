#!/bin/bash
clade=8-genomes
profile=GTDB/$clade/profile.collapse.txt
#profile=GTDB/$clade/profile.txt
tree=GTDB/$clade/bac120.nwk
for alpha in 0.1 0.5 1 2;do
  for lrp in 0.5 1 2;do
  scores=GTDB/$clade/scores.${alpha}.${lrp}.pruned.txt
  sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=8 --nodes=1 --wrap "scripts/scoring-context-conservation.py -p $profile -t $tree --scores $scores -td GTDB/$clade/tmp/${alpha}.${lrp}.pruned --min-tips 2 -lrp $lrp --jobs 16 --prune-tree  > GTDB/$clade/scores.${alpha}.${lrp}.log 2>&1"
  done
done
