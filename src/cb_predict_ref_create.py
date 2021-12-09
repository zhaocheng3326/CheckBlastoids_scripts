# conda activate cb
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
import sys
import anndata
import logging

np.set_printoptions(threshold=200)
pd.set_option("max_rows", 20)
tf.logging.set_verbosity(0)
cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0


wk_dir = sys.argv[1]

ref_pj = ['SPH2016','D3post','CS7']
cb_list={}
if os.path.exists(wk_dir):
    os.chdir(wk_dir)
    for pj in ref_pj :
        cb_list[pj]=cb.data.ExprDataSet.from_anndata(sc.read_h5ad("ref_"+pj+".h5ad"))
        cb_list[pj].uns={"seurat_genes":np.array(cb_list[pj].var.index[cb_list[pj].var['vst.variable']==1])}
        models = []
        for i in range(4):
            models.append(cb.directi.fit_DIRECTi(cb_list[pj], genes=cb_list[pj].uns["seurat_genes"],latent_dim=10, cat_dim=None, random_seed=i))
            blast = cb.blast.BLAST(models, cb_list[pj])
            blast.save(pj+".blast")




