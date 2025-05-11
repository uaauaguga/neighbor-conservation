#!/usr/bin/env python
from scipy.stats import kstest
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

def main():
    idx2pair_id = {}
    with open("GTDB/Gammaproteobacteria.downstream.pairs.txt") as f:
        for line in f:
            idx, pair_id = line.strip().split("\t")
            idx2pair_id[idx] = pair_id
    indir = "GTDB/Gammaproteobacteria.downstream.CTMC"
    scores_by_pair = {}
    background_scores = []
    for txt in tqdm(os.listdir(indir)):
        idx = txt.split(".")[0]
        pair_id = idx2pair_id[idx]
        path = os.path.join(indir, txt)
        table = pd.read_csv(path, sep="\t")
        scores = table["r_gain"]/table["r_loss"]
        background_scores.append(scores.values[-1])
        scores_by_pair[pair_id] = scores.values[0:1000:100]
     
    fout = open("GTDB/Gammaproteobacteria.downstream.scores.txt","w")
    print("pair", "score", "pvalue", sep="\t", file=fout)
    for pair_id in scores_by_pair:
        res = kstest(scores_by_pair[pair_id], background_scores)
        statistic, pvalue = res.statistic, res.pvalue
        print(pair_id, np.mean(scores_by_pair[pair_id]), pvalue, sep="\t",file=fout)
    fout.close()

if __name__ == "__main__":
    main()
