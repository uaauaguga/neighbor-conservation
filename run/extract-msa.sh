#!/bin/bash
for dataset in $(ls GTDB);do
  sem -j 4 python scripts/extract-msa.py -d $dataset > GTDB/${clade}/bac120.tree.log 2>&1
done
sem --wait
