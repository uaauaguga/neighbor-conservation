#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description='summarize context pfams')
parser.add_argument('--input', '-i', type=str, required=True, help='Input context annotations')
parser.add_argument('--output','-o',type=str, required=True, help='Output context pfams')
args = parser.parse_args()

fout = open(args.output,"w")

direction = ["distinct","identical"]
with open(args.input) as f:
    for line in f:
        fields = line[:-1].split("\t")
        #OG0000606       NP_414543.1     NP_414542.1     OG0002928       1       81              NP_414544.1     OG0000607       1       1       PF00288.31,PF08544.18
        group_id, gene_id =  fields[:2]
        upstream_gene_id, downstream_gene_id = fields[2], fields[7]        
        upstream_pfams = fields[6].split(",")
        downstream_pfams = fields[11].split(",")
        gene_1, gene_2 = sorted([upstream_gene_id, gene_id])
        for pfam_id in upstream_pfams:
            if len(pfam_id) > 0:
                print(group_id, gene_id, pfam_id, "upstream", direction[int(fields[4])], sep="\t", file=fout)
        gene_1, gene_2 = sorted([downstream_gene_id, gene_id])
        for pfam_id in downstream_pfams:
            if len(pfam_id) > 0:
                print(group_id, gene_id, pfam_id, "downstream", direction[int(fields[9])], sep="\t", file=fout)
fout.close()
               
