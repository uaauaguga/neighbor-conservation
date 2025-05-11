#!/bin/bash
for bed in $(ls genome/bed );do
  echo $bed
  scripts/gather-contexts.py -i genome/bed/$bed -o output/single-copy-ortholog-gene-contexts/${bed%.*}.txt -a genome/proteins.pfam/${bed%.*}.txt
done
