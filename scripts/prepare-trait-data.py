#!/usr/bin/env python
from collections import defaultdict
import numpy as np
import os
import pandas as pd
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("nex")


def main():    
    parser = argparse.ArgumentParser(description='prepare trait data in nex format')
    parser.add_argument('--input', '-i', required = True, help="input phylogenetic profile")
    parser.add_argument('--output','-o', required = True, help="trait data in nex format")
    parser.add_argument('--names','-n', required = True, help="names of the traits")
    parser.add_argument('--max-nchar','-mn', default=3000, help="max number of characters")
    parser.add_argument('--location','-l', default = "both", 
                        choices = ["upstream","downstream","both"], help="consider pfam in specific location")
    parser.add_argument('--direction','-d', default = "both",
                        choices = ["identical","distinct","both"], help="consider pfam in specific direction")
    args = parser.parse_args()
    traits = defaultdict(list)
    with open(args.input) as f:
        header = next(f).strip().split("\t")
        genome_ids = header[1:] 
        for line in f:
            fields = line.strip().split("\t")
            name = fields[0]
            group_id, location, pfam_id, direction = name.split("-")
            if (args.location != "both") and (location != args.location):
                continue
            if (args.direction != "both") and (direction != args.direction):
                continue
            trait = np.array(fields[1:])
            trait[(trait != '0') & (trait != '1') ] = "?"
            trait = "".join(trait)
            traits[group_id].append((pfam_id, trait))        
    prescences = []
    pairs = []
    encountered_traits = {}
    for group_id in traits:
        traits_by_pfam_ids = defaultdict(list)
        for pfam_id, trait in traits[group_id]:
            traits_by_pfam_ids[pfam_id].append(trait)        
        combined_traits = {}
        for pfam_id in traits_by_pfam_ids:
            combined_trait = np.array(list(traits_by_pfam_ids[pfam_id][0]))
            if len(traits_by_pfam_ids[pfam_id])>0:
                for trait in traits_by_pfam_ids[pfam_id][1:]:
                    trait = np.array(list(trait))
                    #update combined trait
                    combined_trait[trait=='1'] = '1'           
            combined_traits[pfam_id] = combined_trait
        n_informative_max = 0
        for pfam_id in combined_traits:
            n_informative = (combined_traits[pfam_id] != "?").sum()
            if n_informative > n_informative_max:
                pfam_id_max = pfam_id
                n_informative_max = n_informative_max
        combined_trait = combined_traits[pfam_id_max]
        pair = group_id + "-" + pfam_id_max
        combined_trait_s = "".join(combined_trait)
        if combined_trait_s not in encountered_traits:
            encountered_traits[combined_trait_s] = [pair]
        else:
            encountered_traits[combined_trait_s].append(pair)
            continue
        prescences.append(combined_trait)
        pairs.append(pair)

    for combined_trait_s in encountered_traits:
        if len(encountered_traits[combined_trait_s]) > 1:
            print("collapse:", combined_trait_s, ",".join(encountered_traits[combined_trait_s]),sep="\t")        

    n_char_all = len(prescences)
    logger.info(f"Number of characters is {n_char_all} .")
    if n_char_all > args.max_nchar:
        logger.info(f"Sampling {args.max_nchar} characters .")
        indices = np.arange(n_char_all)
        np.random.shuffle(indices)
        indices = indices[:args.max_nchar]
        prescences = [prescences[i] for i in indices]
        pairs = [pairs[i] for i in indices]

    logger.info("Save names of traits finally used .")
    with open(args.names,"w") as f:
        for i, pair in enumerate(pairs):
            print(i+1,pair, file=f, sep="\t")  
          
    prescences = np.array(prescences).T
    fout = open(args.output,"w")
    nchar = prescences.shape[1]
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
    main()
