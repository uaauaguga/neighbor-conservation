#!/usr/bin/env python
import os
import pandas as pd
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("prepare profile")

def main():
    parser = argparse.ArgumentParser(description='prepare phylogenetic profile')
    parser.add_argument('--input-directory','-id', required=True, help="genes with upstream and downstream pfam annotation")
    parser.add_argument('--copy-number', '-cn', required=True, help="copy number of each OG")
    parser.add_argument('--profile', '-p', required = True, help="output phylogenetic profile")
    parser.add_argument('--minimal-informative-fraction', '-mif', default = 0.2, type = float, help="minimal informative fraction required.")
    args = parser.parse_args()


    # load pfam annotations of up/downstream protein flanking 
    # ortholog groups that are mostly single copy genes
    logger.info("Load gene - Pfam pairs ...")
    records = set()
    
    for txt in os.listdir(args.input_directory):
        if not txt.endswith(".txt"):
            continue
        path = os.path.join(args.input_directory, txt)
        genome_id = txt[:-4]
        with open(path) as f:
            for line in f:
                fields = line.strip().split("\t")
                group_id = fields[0]
                pfam_id = fields[2]
                location = fields[3]
                direction = fields[4] 
                pair = "-".join([group_id, location, pfam_id, direction])                
                genome_id = txt[:-4]
                if (pair, genome_id) in records:
                    continue
                records.add((pair, genome_id, 1))

    # gather phylogenetic profile
    records = list(records)
    profile = pd.DataFrame.from_records(records)
    profile.columns = ["pair","genome","prescence"]
    profile = profile.pivot(index="pair",columns="genome",values="prescence").fillna(0).astype(int)
    
    logger.info(f"Load gene copy numbers ...")
    # copy number of genes belong to ortholog group in each genome
    counts = pd.read_csv(args.copy_number,sep="\t",index_col=0)

    # distinguish absent of flanking pfam from absent of otholog gene
    values = counts.loc[profile.index.map(lambda x:x.split("-")[0]),profile.columns].values
    profile.values[values > 1 ] = 2
    profile.values[values == 0 ] = -1
    mask = (profile == 0 ) | (profile == 1 )
    logger.info(f"{profile.shape[0]} genes before filter .")
    cutoff = int(mask.shape[1]*args.minimal_informative_fraction) 
    logger.info(f"The cutoff is set to {cutoff}/{mask.shape[1]} .")
    profile = profile[mask.sum(axis=1) >= cutoff]     
    logger.info(f"{profile.shape[0]} genes passed filter .")
    mask = (profile == 0 ) | (profile == 1 )
    informative_fraction = mask.sum(axis=1).mean()/mask.shape[1]
    logger.info(f"{round(informative_fraction,3)} of the entries on average are informative.")
    positive_fraction = (profile[mask] == 1).sum().sum() / mask.sum().sum()
    logger.info(f"{round(positive_fraction,3)} of the entries on average are positive.")
    logger.info("Saving phylogenetic profile ...")
    profile.to_csv(args.profile,sep="\t")

    logger.info("All done.")

               
if __name__ == "__main__":
    main()
