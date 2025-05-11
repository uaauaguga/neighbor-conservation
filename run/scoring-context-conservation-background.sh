#!/bin/bash
for clade in $(ls GTDB);do
  profile=GTDB/$clade/profile.collapse.txt
  tree=GTDB/$clade/bac120.nwk
  scores=GTDB/$clade/single.emperical.background.pruned.txt
  #[ -s GTDB/$clade/single.emperical.background.pruned.txt ] || 
  sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=4 --nodes=1 --wrap "scripts/scoring-context-conservation.py -p $profile -t $tree --scores $scores -td GTDB/$clade/tmp --prune-tree > GTDB/$clade/single.emperical.background.pruned.log 2>&1"
   scores=GTDB/$clade/single.emperical.background.unpruned.txt
   sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=4 --nodes=1 --wrap "scripts/scoring-context-conservation.py -p $profile -t $tree --scores $scores -td GTDB/$clade/tmp > GTDB/$clade/single.emperical.background.unpruned.log 2>&1"
done
