#!/usr/bin/env python
from collections import defaultdict
import numpy as np
import argparse
from multiprocessing import Pool
import subprocess
import pandas as pd
from scipy.stats import kstest,ks_2samp, ranksums
from arviz import ess
from ete3 import Tree
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('context conservation')

def prepare_nexus(prescence):    
    lines = "" 
    lines += "#NEXUS\n"
    lines += "BEGIN DATA;\n"
    lines += f"DIMENSIONS  NCHAR=1 NTAX={prescence.shape[0]};\n"
    lines += 'FORMAT DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="01";\n'
    lines += "\n"
    lines += "\t MATRIX\n"
    #print(''.join(prescence.astype(str).values))
    for i in range(prescence.shape[0]):
        if prescence.values[i] in [0,1]:
            value = str(prescence.values[i])
        else:
            value = '?'
        line = " ".join(["\t", prescence.index[i], value]) + "\n"
        lines += line
    lines += ";\n"
    lines += "END;\n"
    return lines

def prune(tree, names):
    nodes = []
    for node in tree.traverse():
        name = node.name
        if len(name) == 0:
            continue
        if name in names:
            nodes.append(node)
    tree.prune(nodes)
    return tree

def get_node_names(tree):
    names = []
    for node in tree.traverse():
        names.append(node.name)
    return names

def scoring(nex, nwk, path, alpha, grp, lrp, positive_number, n):
    logger.info(f"run MCMC for {nex} ...")
    cmd = ["rb", "scripts/estimate-rate.free.binary.single.site.Rev", "--args", nex, nwk, path,  str(alpha), str(grp), str(lrp)] 
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    #print(" ".join(cmd))
    table = pd.read_csv(path,sep="\t",index_col=0).reset_index().iloc[-1000:,:]
    r_gains = table["r_gain"].values
    r_losses = table["r_loss"].values
    N = len(r_gains)    
    ESS = int(min(ess(r_gains), ess(r_losses)))
    idx = np.round(np.linspace(0, N - 1, ESS)).astype(int)
    r_gains = r_gains[idx]
    r_losses = r_losses[idx]
    median_r_gain = np.median(r_gains)
    median_r_loss = np.median(r_losses)
    ratio_of_median = median_r_gain/median_r_loss
    ratios = r_gains/r_losses
    median_of_ratio = np.median(ratios)
    character = nex.split("/")[-1][:-4]
    bg_ratio = args.gain_rate_prior/args.loss_rate_prior
    ratio_of_median = ratio_of_median/bg_ratio
    median_of_ratio = median_of_ratio/bg_ratio
    if args.background is not None:
        pvalue1 = (background > ratio_of_median)/ratio_of_median.shape[0]
        pvalue2 = (background > median_of_ratio)/ratio_of_median.shape[0]
        print(character, ratio_of_median, median_of_ratio, positive_number, n,  pvalue1, pvalue2, sep="\t", file=fout)
    else:
        print(character, ratio_of_median, median_of_ratio, positive_number, n, sep="\t", file=fout)
    fout.flush()
    os.remove(nex)
    os.remove(nwk)
    os.remove(path)
    return 0

def main():
    parser = argparse.ArgumentParser(description='estimate turnover rate of surrounding genes')
    parser.add_argument('--profile', '-p', type=str, required=True, help='The phylogenetic profiles')
    parser.add_argument('--tree', '-t', type=str, required=True, help='The phylogenetic tree')
    parser.add_argument('--alpha-prior', '-ap', type=float, default=2, help='prior for shape of gamma distribution')
    parser.add_argument('--gain-rate-prior', '-grp', type=float, default=0.2, help='gain rate prior')
    parser.add_argument('--loss-rate-prior', '-lrp', type=float, default=2.0, help='loss rate prior')
    parser.add_argument('--scores', '-s', type=str, required=True, help='where to save the scores')
    parser.add_argument('--temp-dir', '-td', type=str, required=True, help='temp directory')
    parser.add_argument('--prune-tree', '-pt', action = "store_true", help='whether exclude unobserved tips')
    parser.add_argument('--min-tips', '-mt', type=int, default=10, help='minimal number of observed tips')
    parser.add_argument('--background', '-bg', type=str, help='background scores')
    parser.add_argument('--jobs', '-j', type=int, default=16, help='number of wokers to use')
    parser.add_argument('--min-positive-fraction', '-mpf', type=float, default=0, help='number of wokers to use')
    global args
    args = parser.parse_args()

    logger.info("Load phylogenetic profile ...")
    profile = pd.read_csv(args.profile,sep="\t",index_col=0)    
    
    logger.info("Load tree ...")
    tree = Tree(newick=args.tree,format=1,quoted_node_names=True) 


    genome_ids = get_node_names(tree)
    common_genome_ids = np.intersect1d(genome_ids,profile.columns)
    


    logger.info(f"{len(genome_ids)} genomes in tree, {profile.shape[1]} in profile,")
    logger.info(f"{len(common_genome_ids)} are common.")
    
    global background
    if args.background is not None:
        background = np.array(open(args.background).read().strip().split("\n")).astype(float)

    tree = prune(tree, common_genome_ids)
    profile = profile.loc[:,common_genome_ids]
    

    logger.info("Run MCMC analysis for rate estimation ...")

    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)

    global fout

         

    fout = open(args.scores,"w")
    occured_profiles = set()   
    pool = Pool(args.jobs) 
    workers = []
    nchar = 0
    for character in profile.index:
        prescence = profile.loc[character]
        mask = (prescence.values == 0) | (prescence.values == 1)
        n, N = mask.sum(), mask.shape[0]
        if n < args.min_tips:
            continue
        p = "".join(prescence.astype(str).values)
        if p in occured_profiles:
            continue
        occured_profiles.add(p)
        used_tree = tree.copy()
        genome_ids = np.array(list(profile.columns))[mask]
        if args.prune_tree:
            used_tree = prune(used_tree, set(list(genome_ids)))
            prescence = prescence.loc[genome_ids]
        positive_number = (prescence==1).sum()
        if positive_number/n < args.min_positive_fraction:
            continue
        print(f"{n} in {N} tips are observed, number of positive is {positive_number}.")

        nwk = os.path.join(args.temp_dir,character + ".nwk")
        with open(nwk,"w") as f:
            f.write(used_tree.write(format=1)) 
        nex = os.path.join(args.temp_dir,character + ".nex")        
        with open(nex,"w") as f:
            f.write(prepare_nexus(prescence))
        path = os.path.join(args.temp_dir,character + ".txt")
        #scoring(nex, nwk, path, args.alpha_prior, args.gain_rate_prior, args.loss_rate_prior, positive_number/N)       
        workers.append(pool.apply_async(func=scoring,args=(nex, nwk, path, args.alpha_prior, args.gain_rate_prior, args.loss_rate_prior, positive_number, n))) 
        nchar += 1
    for worker in workers:
        _ = worker.get()
    fout.close()
                        



if __name__ == "__main__":
    main()
