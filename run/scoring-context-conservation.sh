#!/bin/bash
#clade=8-genomes
clade=Acidobacteriae
profile=GTDB/$clade/profile.collapse.txt
#profile=GTDB/$clade/profile.txt
tree=GTDB/$clade/bac120.nwk
scores=GTDB/$clade/scores.txt
scripts/scoring-context-conservation.py -p $profile -t $tree --scores $scores -td GTDB/$clade/tmp --min-tips 2 
