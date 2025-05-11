#!/usr/bin/env python
import pandas as pd
from random import shuffle
import os
from pyfaidx import Fasta
from Bio.Seq import Seq

def get_taxonomy(s, l = "c"):
    t = "."
    for c in s.split(";"):
        if c.startswith(l):
            t = c[3:] 
    return t

def main():
    #gtdb_taxonomy
    level = "f"
    #level = "c"
    table = pd.read_csv("../sRNA-annotation/bac120_metadata_r207.rep.tsv",sep="\t",index_col=0)
    table[level] = table["gtdb_taxonomy"].map(lambda x:get_taxonomy(x, level))
    sizes = table.groupby(level).apply(lambda x:x.shape[0]) 
    sizes = sizes.sort_values(ascending=False)
    tax_ids = sizes.index[:20]
    
    sampled_genome_ids_by_clade = {}
    for tax_id in tax_ids:
        genome_ids = list(table[table[level] == tax_id].index.map(lambda x:x[x.find("_")+1:]))
        shuffle(genome_ids)
        sampled_genome_ids_by_clade[tax_id] = genome_ids[:64]
    
    genome_id2chunk_id = {}
    for txt in os.listdir("/data2/lulab1/jinyunfan/sRNA-annotation/data/GTDB/genome-ids"):
        chunk_id = txt[:-4]
        genome_ids = open("/data2/lulab1/jinyunfan/sRNA-annotation/data/GTDB/genome-ids/" + txt).read().strip().split("\n")
        for genome_id in genome_ids:
            genome_id2chunk_id[genome_id] = chunk_id

    
    for tax_id in sampled_genome_ids_by_clade:
        print(f"processing {tax_id} ...")
        if not os.path.exists(f"GTDB/{tax_id}"):
            os.mkdir(f"GTDB/{tax_id}")
            os.mkdir(f"GTDB/{tax_id}/proteins")
            os.mkdir(f"GTDB/{tax_id}/genes")
        genome_ids = sampled_genome_ids_by_clade[tax_id]
        for genome_id in genome_ids:
            chunk_id = genome_id2chunk_id[genome_id]
            fasta = Fasta(f"/data2/lulab1/jinyunfan/sRNA-annotation/data/GTDB/fasta/{chunk_id}/{genome_id}.fa")
            bed =  f"/data2/lulab1/jinyunfan/sRNA-annotation/data/GTDB/bed/{genome_id}.bed"
            fgene = open(f"GTDB/{tax_id}/genes/{genome_id}.bed","w")
            fprotein = open(f"GTDB/{tax_id}/proteins/{genome_id}.faa","w")
            with open(bed) as f:
                for line in f:
                    fields = line.strip().split("\t")
                    seq_id, start, end = fields[:3]
                    protein_id = fields[3] 
                    strand = fields[5] 
                    start, end = int(start), int(end)
                    sequence = fasta[seq_id][start:end]
                    if strand == "-":
                        sequence = sequence.reverse.complement
                    sequence = str(sequence)
                    aas = Seq(sequence).translate()                 
                    if aas.endswith("*"):
                        aas = aas[:-1]
                    fgene.write(line)
                    print(f">{protein_id}",file=fprotein)
                    print(aas,file=fprotein)
            fprotein.close()     
            fgene.close()
        
if __name__ == "__main__":
    main()
