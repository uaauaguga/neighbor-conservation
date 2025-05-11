#!/usr/bin/env python
import argparse
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("prepare profile")
import subprocess

def main():
    parser = argparse.ArgumentParser(description='extract msa')
    parser.add_argument('--dataset', '-d', required=True, type = str, help="dataset to use")
    args = parser.parse_args()
    clade = args.dataset 
    genome_ids = set()
    logger.info("count number of genomes ...")
    for bed in os.listdir(f"GTDB/{clade}/genes"):
        genome_ids.add(bed[:-4])
    
    fout = open(f"GTDB/{clade}/bac120.faa","w")

    logger.info(f"extract from msa ...") 
    n = 0
    with open("gtdbtk.bac120.msa.fasta") as f:
        for header in f:
            header2 = ">" + header[4:]
            sequence = next(f)
            genome_id = header2[1:].split(" ")[0]
            if genome_id in genome_ids:
                n += 1
                print(header2.strip(), file=fout)
                print(sequence.strip(), file=fout)
    fout.close()

    #stdout = subprocess.DEVNULL,
    #stderr = subprocess.STDOUT 
    fout = open(f"GTDB/{clade}/bac120.nwk","w")
    logger.info(f"see {n} in {len(genome_ids)} sequences.")
    logger.info("Build tree ...")
    subprocess.run(["FastTree", "-wag", '-gamma', f"GTDB/{clade}/bac120.faa"], stdout = fout)
    fout.close()
    logger.info("All done.")
                
if __name__ == "__main__":
    main()
