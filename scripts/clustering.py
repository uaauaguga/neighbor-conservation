#!/usr/bin/env python
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances


def main():
    profile = pd.read_csv("GTDB/Clostridia/profile.0.4.txt",sep="\t",index_col=0)
    ac = AgglomerativeClustering(distance_threshold=0.01,n_clusters=None) 
    ac.fit(profile.values)
    print(ac.labels_.shape) 
    print(ac.labels_[:50])
    


if __name__ == "__main__":
    main()
