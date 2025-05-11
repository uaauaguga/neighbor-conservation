#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='get closest gene pairs')
    parser.add_argument('--input', '-i', type=str, required=True, help='input table')
    parser.add_argument('--output', '-o', type=str, required=True, help='output table')
    args = parser.parse_args()
    #NC_000913.3     190     255     NP_414542.1     RDBECOLIOPC00853        +       b0001   thrL
    last_seq_id = ""
    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end, protein_id, operon_id, strand = fields[:6]
            if (protein_id == ".") or (operon_id == "."):
                continue
            if seq_id == last_seq_id:
                same_operon = int(last_operon_id == operon_id)
                distance = int(start) - int(last_end)
                same_strand = int(last_strand == strand)
                print(last_protein_id, protein_id, same_operon, same_strand, distance, sep="\t", file=fout)
            last_seq_id, last_start, last_end, last_protein_id, last_operon_id, last_strand = seq_id, start, end, protein_id, operon_id, strand
    fout.close()
if __name__ == "__main__":
    main()
