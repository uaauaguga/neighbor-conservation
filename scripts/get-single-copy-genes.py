#!/usr/bin/env python
import pandas as pd
def main():
    parser = argparse.ArgumentParser(description='summarize single copy orthologs')
    parser.add_argument('--input-directory', '-id', required=True, help="input query proteins")
    parser.add_argument('--output','-o', required=True, help="output hits")
    args = parser.parse_args()
    table = pd.read_csv("{args.input_directory}/Orthogroups/Orthogroups.GeneCount.tsv",sep="\t",index_col=0)
    table = table[(table > 1).iloc[:,:-1].sum(axis=1) == 0]
    group_ids = table = table[table.iloc[:,:-1].sum(axis=1) >=5].index
    table = pd.read_csv("genome/orthologs/Results_Aug09/Orthogroups/Orthogroups.tsv",sep="\t",index_col=0)
    table = table.loc[group_ids,:]
    table.to_csv("genome/orthologs.txt",sep="\t") 
    
if __name__ == "__main__":
    main()
