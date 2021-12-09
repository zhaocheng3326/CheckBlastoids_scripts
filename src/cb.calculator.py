#!/usr/bin/env python


"""
calculate cell blast based on SPH2016, 3Dpost, CS7 dataset
Written by Cheng Zhao 2021
"""
import os, sys, re

from argparse import ArgumentParser
parser = ArgumentParser(description='calculate cb from h5ad file')

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise string+" ,dir not exists"


# [Required input]
parser.add_argument('-i', '--inFile', metavar='inFile', default=[], help="Input h5ad file", required=True)
parser.add_argument('-c', '--cbDatabase', metavar='cbDatabase', default=[], help="dir for cbDatabase", required=True,type=dir_path)
parser.add_argument('-o', '--cbOutput', metavar='cbOutput', default=[], help="output file", required=True)

args = parser.parse_args()

# imprt libarary
import time
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import pandas as pd
import tensorflow as tf
import scanpy as sc
import Cell_BLAST as cb
import os
import anndata
import logging


cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0

input_file = args.inFile
cb_dir = args.cbDatabase

#output_file= input_file +".query.ref.csv"
output_file= args.cbOutput
# loading input file
query_file=cb.data.ExprDataSet.from_anndata(sc.read_h5ad(input_file))
query_file.uns={"seurat_genes":np.array(query_file.var.index[query_file.var['vst.variable']==1])}


query_hits=[]
for ref in ['SPH2016','D3post','CS7']:
    blast = cb.blast.BLAST.load(cb_dir+"/"+ref+".blast")
    temp=blast.query(query_file).reconcile_models().filter(by="pval", cutoff=0.1).to_data_frames()
    for x in temp:
        temp[x]['ref_cell']=temp[x].index
        temp[x]['ref_cell']=temp[x].index
        temp[x]['query_cell']=x
        temp[x]=temp[x][['ref_cell','query_cell','hits','dist','pval']]
        query_hits.append(temp[x])
    print("[INFO] Done with "+ref+" prediction")
pd.concat(query_hits).to_csv(output_file,sep="\t",index=False)
print("[INFO] Done with cell_blast prediction")
