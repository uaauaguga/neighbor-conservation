#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("homolog search")
import subprocess
import os
from collections import defaultdict
from shutil import rmtree

def main():
    parser = argparse.ArgumentParser(description='protein homolog search')
    parser.add_argument('--query', '-q', required=True, help="input query proteins")
    parser.add_argument('--output','-o', required=True, help="output hits")
    parser.add_argument('--target','-t', default="/data2/lulab1/jinyunfan/Pfam/pfam_profile", help="target database")
    parser.add_argument('--tempdir','-td', required=True, help="temporary directory")
    args = parser.parse_args()

    if not os.path.exists(args.tempdir):
        os.mkdir(args.tempdir)

    tdb = args.target

    mmseqs = "/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin/mmseqs" 
    qdb = os.path.join(args.tempdir, "query")
    if not os.path.exists(qdb + ".dbtype"):
        logger.info("build protein mmseqs database ...")
        cmd = [mmseqs, "createdb", "--dbtype", "1", args.query, qdb ]
        subprocess.run(cmd)
    else:
        logger.info("protein db exists .")

    hdb = os.path.join(args.tempdir, "hits")
    if not os.path.exists(hdb+".dbtype"):
        logger.info("homolog search ...")
        cmd = [mmseqs, "search", "-s", "2", "--max-seqs", "500", qdb, tdb, hdb, args.tempdir]
        subprocess.run(cmd)
    else:
        logger.info("hit db exists .")


    tsv = os.path.join(args.tempdir, "hits.tsv")
    if not os.path.exists(tsv):
        logger.info("format homolog search hits ...")
        cmd = [mmseqs, "convertalis", "--format-mode", "2", qdb, tdb, hdb, tsv ]
        subprocess.run(cmd)
    else:
        logger.info("hit file exists .")


    logger.info("fltering hits ...")
    fout = open(args.output,"w")
    with open(tsv) as f:
        for line in f:
            fields = line.strip().split("\t")    
            aligned_length = int(fields[9]) - int(fields[8])
            coverage = aligned_length/int(fields[13])                   
            fout.write(line)
    fout.close()

    logger.info("remove temp files ...")
    rmtree(args.tempdir)    

    logger.info("all done.")

    
    

              

if __name__ == "__main__":
    main()
