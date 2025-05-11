#!/bin/bash
indir=output/single-copy-ortholog-gene-contexts
outdir=output/single-copy-ortholog-gene-contexts-pfam
for txt in $(ls $indir);do
  scripts/summarize-context-pfam.py -i $indir/$txt -o $outdir/$txt
done
