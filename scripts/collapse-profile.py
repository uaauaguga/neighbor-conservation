#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='collapse profile')
    parser.add_argument('--input','-i', type=str , required = True, help="Input profile")
    parser.add_argument('--output','-o', type=str, required=True, help="Output profile")
    args = parser.parse_args()
    table = pd.read_csv(args.input,sep="\t",index_col=0)

    def extract_group(x):
        OG_id, location, pfam_id, direction = x.split("-")
        return OG_id + "-" + pfam_id

    table["group"] = table.index.map(extract_group)

    def combine(x):
        names = x.columns[:-1]
        x = x.values[:,:-1]
        if x.shape[0] == 1:
            d = x[0,:]
        else:
            d = np.zeros(x.shape[1])
            d[(x==2).sum(axis=0)>0] = 2
            d[(x==-1).sum(axis=0)>0] = -1
            d[(x==1).sum(axis=0)>0] = 1
        return pd.Series(data=d,index=names).astype(int)

    combined_table = table.groupby("group").apply(combine)
    combined_table.to_csv(args.output,sep="\t")
if __name__ == "__main__":
    main()
