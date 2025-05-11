#!/usr/bin/env python
import pandas as pd
import numpy as np
def main():
    clade_ids = open("clades.txt").read().strip().split("\n")
    scores = []
    for clade_id in clade_ids:
        table = pd.read_csv(f"GTDB/{clade_id}/single.emperical.background.pruned.txt",sep="\t",header=None)
        scores += list(table[1].values)
    scores = np.array(scores)
    #for i in np.linspace(0,20,2001):
    #   print(i,(scores>i).sum()/scores.shape[0]) 
    np.random.shuffle(scores)
    with open("emperical-background.txt","w") as f:
        for score in scores:
            print(score, file=f)     
   
if __name__ == "__main__":
    main()
