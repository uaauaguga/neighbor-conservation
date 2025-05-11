#!/usr/bin/env python
from Bio import Entrez
def fetch_url(seq_id):
    Entrez.email = 'jyf20@mails.tsinghua.edu.cn'
    handle = Entrez.esearch(db="assembly", term=seq_id, retmax='1')
    record = Entrez.read(handle)
    if len(record['IdList']) == 0:
        return ".", "."
    esummary_handle = Entrez.esummary(db="assembly", id=record['IdList'][0], report="full")
    esummary_record = Entrez.read(esummary_handle)
    url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    organism = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
    return url, organism

if __name__ == "__main__":
    seq_ids = open("seq-ids.txt").read().strip().split("\n")
    for seq_id in seq_ids:
        url, organism = fetch_url(seq_id)
        print(seq_id, url, organism, sep="\t")
       

