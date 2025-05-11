#!/usr/bin/env python
import argparse
import pandas as pd
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("annotate context")

def main():
    parser = argparse.ArgumentParser(description='gather contexts of selected proteins')
    parser.add_argument('--input', '-i', type=str, required=True, help='input CDS in bed format')
    parser.add_argument('--output', '-o', type=str, required=True, help='context of focused genes')
    parser.add_argument('--pfam-annotation', '-a', type=str, required=True, help='pfam annotations')
    parser.add_argument('--og-annotation', '-g', type=str, required=True, help='ortholog annotation')
    args = parser.parse_args()

    logger.info("Load single copy genes and corresponds ortholog groups ...")
    genome_id = args.input.split("/")[-1][:-4]

    genes = pd.read_csv(args.og_annotation,sep="\t",index_col=0)
    genes = genes[genome_id]
    genes = genes[~genes.isna()]
    genes = genes[~genes.map(lambda x:"," in x)]

    gene_id2group_id = {}
    for group_id, gene_id in zip(genes.index, genes.values):
        gene_id2group_id[gene_id] = group_id


    logger.info("Load gene order ...")
    records = []    
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end, protein_id = fields[:4]
            start, end = int(start), int(end)
            strand = fields[-1]
            records.append((seq_id, start, end, protein_id, strand))

    logger.info("Load pfam annotation ...")
    pfam_hits = defaultdict(set)
    with open(args.pfam_annotation) as f:
        for line in f:
            fields = line.strip().split("\t")
            protein_id, pfam_id = fields[:2]
            pfam_hits[protein_id].add(pfam_id)

    logger.info("Saving annotated single copy genes ...")
    fout = open(args.output,"w")    
    for i in range(1,len(records)-1):
        seq_id, start, end, protein_id, strand = records[i]
        if protein_id in gene_id2group_id:
            group_id = gene_id2group_id[protein_id]
            last_seq_id, last_start, last_end, last_protein_id, last_strand = records[i-1]
            next_seq_id, next_start, next_end, next_protein_id, next_strand = records[i+1]        
            upstream_distance = start - last_end
            upstream_direction = int(last_strand == strand)    
            downstream_distance = next_start - end
            downstream_direction = int(strand == next_strand)
            if strand == "-":
                last_seq_id, last_start, last_end, last_protein_id, last_strand = records[i+1]
                next_seq_id, next_start, next_end, next_protein_id, next_strand = records[i-1]
                upstream_distance, downstream_distance = downstream_distance, upstream_distance
                upstream_direction, downstream_direction = downstream_direction, upstream_direction 
            print(group_id, protein_id, 
                  last_protein_id, gene_id2group_id.get(last_protein_id,"."), upstream_direction, upstream_distance, ",".join(sorted(list(pfam_hits[last_protein_id]))),
                  next_protein_id, gene_id2group_id.get(next_protein_id,"."),downstream_direction, downstream_distance, ",".join(sorted(list(pfam_hits[next_protein_id]))), file=fout, sep="\t")
    fout.close()
if __name__ == "__main__":
    main()
