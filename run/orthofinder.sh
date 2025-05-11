#!/bin/bash
for clade in $(ls GTDB);do
  sbatch --exclude=node[15-25] -p Acluster --ntasks-per-node=10 --nodes=1 --wrap="orthofinder -t 16 -f GTDB/$clade/proteins -o GTDB/$clade/orthologs > GTDB/$clade/orthologs.log 2>&1"
done
