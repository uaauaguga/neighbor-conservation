#!/usr/bin/env python
import pandas as pd
def main():
    table = pd.read_csv("genome/orthologs/Results_Aug09/Orthogroups/Orthogroups.GeneCount.tsv",sep="\t",index_col=0)
    table = table[(table > 1).iloc[:,:-1].sum(axis=1) == 0]
    table = table[table["GCF_000005845.2_ASM584v2_protein"]>0]
    group_ids = table = table[table.iloc[:,:-1].sum(axis=1) >=5].index
    #int(len(group_ids))
    table = pd.read_csv("genome/orthologs/Results_Aug09/Orthogroups/Orthogroups.tsv",sep="\t",index_col=0)
    table = table.loc[group_ids,:]
    table.to_csv("genome/orthologs.txt",sep="\t") 
    
if __name__ == "__main__":
    main()
