#!/bin/bash
for clade in $(ls GTDB);do
   [ -s GTDB/$clade/profile.collapse.txt ] || sbatch --exclude='node[25-36]' -p Acluster --ntasks-per-node=5 --nodes=1 --wrap "scripts/collapse-profile.py -i GTDB/$clade/profile.txt -o  GTDB/$clade/profile.collapse.txt"
done
