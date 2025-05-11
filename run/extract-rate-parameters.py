#!/usr/bin/env python
from glob import glob
import os
import pandas as pd

def main():
    fout = open("rate.under.differen.context.txt","w")
    print("clade","alpha","gain","loss","direction","location",sep="\t",file=fout)
    for clade in os.listdir("GTDB"):
        for d in ["identical","distinct","both"]:
            for l in ["upstream","downstream","both"]:
                path = f"GTDB/{clade}/rate.{l}.{d}.txt"
                if not os.path.exists(path):
                    continue
                table = pd.read_csv(path, sep="\t", index_col=0)
                #alpha   r_gain  r_loss
                alpha = table["alpha"].iloc[-1000:,].mean()
                r_gain = table["r_gain"].iloc[-1000:,].mean()
                r_loss = table["r_loss"].iloc[-1000:,].mean()
                print(clade,alpha, r_gain,r_loss,d,l,sep="\t",file=fout)
    fout.close()

if __name__ == "__main__":
    main()
