#!/usr/bin/env python
from collections import defaultdict
import numpy as np
def prepare_trait(position,clade):    
    traits = defaultdict(set)
    prescences = []   
    #clade = "Gammaproteobacteria"
    #clade = "Campylobacteria" 
    #clade = "8-genomes"
    clade = "Cyanobacteriia"
    pair_ids = []
    with open(f"GTDB/{clade}.{position}-profile.txt") as f:
        header = next(f).strip().split("\t")
        genome_ids = header[1:] 
        for line in f:
            fields = line.strip().split("\t")
            group_id = fields[0].split("-")[0]            
            if int(fields[-1]) > -10 :
                name = fields[0]
                trait = np.array(fields[1:])
                trait[(trait != '0') & (trait != '1') ] = "?"
                trait = "".join(trait) 
                if trait in traits[group_id]:
                    continue
                if trait.count("?")/len(trait) > 0.1:
                    continue
                traits[group_id].add(trait)
                prescences.append(list(trait))
                pair_ids.append(name)

    with open(f"GTDB/{clade}.{position}.pairs.txt","w") as ft:
        for i,pair_id in enumerate(pair_ids):
            print(i+1, pair_id, sep="\t", file=ft)
    prescences = np.array(prescences).T
    fout = open(f"GTDB/{clade}.{position}.nex","w")
    nchar = min(5000,prescences.shape[1])
    print("#NEXUS", file=fout)
    print("BEGIN DATA;", file=fout)
    print(f"DIMENSIONS  NCHAR={nchar} NTAX={prescences.shape[0]};", file=fout)
    print('FORMAT DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="01";', file=fout)
    print("", file=fout)
    print("\t","MATRIX", file=fout)
    for i in range(prescences.shape[0]):
        print("\t", genome_ids[i], "".join(prescences[i])[:nchar], file=fout)
    print("\n;", file=fout)
    print("END;", file=fout)
    fout.close()

if __name__ == "__main__":
    prepare_trait("upstream","Gammaproteobacteria")
    prepare_trait("downstream","Gammaproteobacteria")
