#!/usr/bin/env python
import pandas as pd
#profile = pd.read_csv("random-profile.txt",sep="\t",index_col=0)
profile = pd.read_csv("random-profile-8-16.txt",sep="\t",index_col=0)
#scores = pd.read_csv("random-profile.scores.txt",sep="\t",index_col=0,header=None).sort_values(by=1,ascending=False)
scores = pd.read_csv("random-profile-8-16.scores.txt",sep="\t",index_col=0,header=None).sort_values(by=1,ascending=False)
p = int(scores.shape[0]/2)
n = 4
profile = profile.loc[list(scores.index[:n]) + list(scores.index[p:p+n]) + list(scores.index[-n:])].T
print(profile.shape)
profile.to_csv("random.profile-8-16.itol.txt")

