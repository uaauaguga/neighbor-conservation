#!/usr/bin/env python
import os
import pandas as pd
def main():
    upstream_pfam_within = {}
    downstream_pfam_within = {}
    with open("output/single-copy-ortholog-gene-contexts-ecoli-operon/context-pfam.txt") as f:
        for line in f:
            fields = line.strip().split("\t")
            #OG0000606       NP_414543.1     PF00288.31      downstream      1
            group_id = fields[0]
            pfam_id = fields[2]
            if fields[3] == "upstream":
                upstream_pfam_within[group_id + "-" + pfam_id] = int(fields[4])
            else:
                downstream_pfam_within[group_id + "-" + pfam_id] = int(fields[4])


    upstream_records = []
    downstream_records = []
    for txt in os.listdir("output/single-copy-ortholog-gene-contexts-pfam"):
        path = os.path.join("output/single-copy-ortholog-gene-contexts-pfam",txt)
        with open(path) as f:
            for line in f:
                fields = line.strip().split("\t")
                group_id = fields[0]
                pfam_id = fields[2]
                direction = fields[4] 
                pair = group_id + "-" + pfam_id
                #OG0000606       NP_414543.1     PF00288.31      downstream
                if fields[3] == "upstream":
                    #if pair in upstream_pfam_within:
                    upstream_records.append((pair + "-" + direction, txt[:-4], 1))
                else:              
                    #if pair in downstream_pfam_within:
                    downstream_records.append((pair + "-" + direction, txt[:-4], 1))

    upstream_profile = pd.DataFrame.from_records(upstream_records)
    upstream_profile.columns = ["pair","genome","prescence"]
    upstream_profile = upstream_profile.pivot(index="pair",columns="genome",values="prescence").fillna(0).astype(int)
    upstream_profile["same operon in E.coli"] = upstream_profile.index.map(lambda x:upstream_pfam_within.get(x[:x.rfind("-")],-1))
    upstream_profile.to_csv("upstream-profile.txt",sep="\t")


    downstream_profile = pd.DataFrame.from_records(downstream_records)
    downstream_profile.columns = ["pair","genome","prescence"]
    downstream_profile = downstream_profile.pivot(index="pair",columns="genome",values="prescence").fillna(0).astype(int)
    downstream_profile["same operon in E.coli"] = downstream_profile.index.map(lambda x:downstream_pfam_within.get(x[:x.rfind("-")],-1))
    downstream_profile.to_csv("downstream-profile.txt",sep="\t")
               
if __name__ == "__main__":
    main()
