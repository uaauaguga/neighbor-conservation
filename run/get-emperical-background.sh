#!/bin/bash
for clade in ls $(GTDB);do
  scripts/scoring-context-conservation.py --profile GTDB/Gammaproteobacteria.downstream-profile.txt --tree GTDB/Gammaproteobacteria/bac120.nwk --scores emperical-background.txt --temp-dir tmp -grp 0.1 -lrp 1.0 -ap 2 -mt 50
done
