# conda cb
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.manifold import TSNE

import matplotlib
import matplotlib.pyplot as plt

# Data from hEP-structures
os.chdir('/home/chenzh/My_project/JP_project/')


Goetz_D = sc.read_10x_mtx("data/predict_data/Sozen_NC_2021/EPSC")

# Data from D5/D6 natural embryos
Goetz_D2 = sc.read_10x_mtx("data/predict_data/Sozen_NC_2021/Embryo")

adata = Goetz_D.concatenate(Goetz_D2, index_unique=None)

adata.var["gene_name"] = adata.var.index
adata.var


# Barcodes/LMOs for EPSCs, D5 hEP-structures, and D6 hEP-structures.  (Multiplexed samples, hence LMO tags)
tags_ep = pd.read_csv("data/predict_data/Sozen_NC_2021/EPSC/GSM5387817_EPSC_LMO-tags.csv")

# Barcodes/LMOs for EPSCs, D5 hEP-structures, and D6 hEP-structures.  (Multiplexed samples, hence LMO tags)
tags_nat = pd.read_csv("data/predict_data/Sozen_NC_2021/Embryo/barcodes.tsv.gz", names=["cell_barcode", "sample_id", "sample_number"], header=None)
tags_nat["sample_id"] = "Nat"
tags_nat["sample_number"] = 7

# Combine the two dataframes, note that order must match order of adata file. (EPSC then natural samples)
tag = (tags_ep, tags_nat)
tags = pd.concat(tag, ignore_index=True)
tags = tags.set_index("cell_barcode").rename_axis(index=None, axis=1)
tags

tags['cell_group'] = [tags['sample_id'][i].split('-')[0] for i in range(len(tags))]
tags.index.name = None

tags


adata.obs = pd.concat([adata.obs, tags], axis=1)

adata.obs


adata = adata[adata.obs['sample_id']!= "unknown"]
adata.obs

# Give labels a more descriptive name
adata.obs["sample_group"] = " "
adata.obs.loc[(adata.obs["cell_group"]=="Nat"), "sample_group"] = "Natural human embryo"
adata.obs.loc[(adata.obs["cell_group"]=="D5"), "sample_group"] = "Day 5 hEP-structures"
adata.obs.loc[(adata.obs["cell_group"]=="D6"), "sample_group"] = "Day 6 hEP-structures"
adata.obs.loc[(adata.obs["cell_group"]=="2D"), "sample_group"] = "hEPSCs in 2D"
# adata.obs.loc[(adata.obs["cell_group"]=="unknown"), "sample_group"] = "unknown"

adata.obs = adata.obs.drop(columns=["batch", "sample_number"])
adata.obs



# Add number of genes and number of counts to adata.obs
adata.obs["n_genes"] = np.asarray((adata.X>0).sum(axis=1)).reshape(-1)
adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).reshape(-1)

# Desginates mitochondrial genes
mito_genes = adata.var.gene_name.str.startswith('MT-')

# For each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)


#examine mitochondrial content
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

matplotlib.rcParams.update({'font.size': 7})
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], groupby="sample_group",jitter=0.4, multi_panel=True)

adata.raw = adata

# Generate list of highly variable genes.
sc.pp.highly_variable_genes(adata,min_mean=0.01,max_mean=8,min_disp=1,n_top_genes=2000,flavor="seurat_v3",n_bins=20)

# Normalize the data.
sc.pp.normalize_total(adata, target_sum=1e4)

# Log the data.
sc.pp.log1p(adata)

# Scale the data.
sc.pp.scale(adata, max_value=10)

# Performs PCA.
sc.tl.pca(adata,svd_solver='arpack', use_highly_variable=True)

# Performs neighbor analysis.
sc.pp.neighbors(adata,n_neighbors=30, n_pcs=50, random_state=42)

sc.tl.louvain(adata)

sc.tl.umap(adata)

# Reads in a list of genes for each lineage from from supplementary table 12 of Liu et al.
df_lineage = pd.read_csv('data/predict_data/Sozen_NC_2021/hEP-structures_MZG/data/supp12.csv')
epi_markers = df_lineage.loc[df_lineage["type"]=="ALL-EPI"]["geneName"].values
hypo_markers = df_lineage.loc[df_lineage["type"]=="ALL-PE"]["geneName"].values
te_markers = df_lineage.loc[df_lineage["type"]=="ALL-TE"]["geneName"].values

# Reads in a list of genes for hEPSCs from Yan et al.
df_epscs = pd.read_csv("data/predict_data/Sozen_NC_2021/hEP-structures_MZG/data/epsc_genes.csv", usecols=[0])
epsc_markers = df_epscs["Genes"].values

# Makes sure that all genes in lineage gene lists are also present in the data set.
epi_markers_present = []
for marker in epi_markers:
    if marker in adata.var["gene_name"]:
        epi_markers_present.append(marker)

hypo_markers_present = []
for marker in hypo_markers:
    if marker in adata.var["gene_name"]:
        hypo_markers_present.append(marker)

te_markers_present = []
for marker in te_markers:
    if marker in adata.var["gene_name"]:
        te_markers_present.append(marker)

epsc_markers_present = []
for marker in epsc_markers:
    if marker in adata.var["gene_name"]:
        epsc_markers_present.append(marker)

# Makes a library for each list of genes.
marker_genes_dict = {"epi":epi_markers_present,
                     "hypo":hypo_markers_present,
                     "te":te_markers_present,
                     "epsc": epsc_markers_present}

# Calculates a "score" for every cell and every lineage, as well as an EPSC score.

ctrl_size_var = 50

sc.tl.score_genes(adata, te_markers_present, ctrl_size=ctrl_size_var, score_name="te_score", use_raw=False)
sc.tl.score_genes(adata, epi_markers_present, ctrl_size=ctrl_size_var, score_name="epi_score", use_raw=False)
sc.tl.score_genes(adata, hypo_markers_present, ctrl_size=ctrl_size_var, score_name="hypo_score", use_raw=False)
sc.tl.score_genes(adata, epsc_markers_present, ctrl_size=ctrl_size_var, score_name="epsc_score", use_raw=False)


# All cells start as undefined.
adata.obs["lineage"] = "undefined"

# Assigns lineage based on where the highest lineage score comes from.
# Assigns EPI group
adata.obs.loc[(adata.obs["epi_score"]>adata.obs["hypo_score"])
              & (adata.obs["epi_score"]>adata.obs["te_score"]), "lineage"] = "epi"

# Assigns HYPO group
adata.obs.loc[(adata.obs["hypo_score"]>adata.obs["te_score"])
              & (adata.obs["hypo_score"]>adata.obs["epi_score"]), "lineage"] = "hypo"

# Assigns TE group
adata.obs.loc[(adata.obs["te_score"]>adata.obs["hypo_score"])
              & (adata.obs["te_score"]>adata.obs["epi_score"]), "lineage"] = "te"


# Shows counts for all lineages.
adata.obs.lineage.value_counts()

# Assigns cells as undefined if lineage is not well defined.
adata.obs.loc[((adata.obs["hypo_score"]<0.08)
               & (adata.obs["te_score"]<0.08)
               & (adata.obs["epi_score"]<0.08)) , "lineage"] ="undefined"
adata.obs.lineage.value_counts()

adata.obs["lineage_id"] = "Undefined"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="Nat")), "lineage_id"] = "Epiblast"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="D5")), "lineage_id"] = "D5 ELCs"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="D6")), "lineage_id"] = "D6 ELCs"

adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="Nat")), "lineage_id"] = "Hypoblast"
adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="D5")), "lineage_id"] = "D5 HLCs"
adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="D6")), "lineage_id"] = "D6 HLCs"

adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="Nat")), "lineage_id"] = "Trophectoderm"
adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="D5")), "lineage_id"] = "D5 TLCs"
adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="D6")), "lineage_id"] = "D6 TLCs"

adata.obs.loc[((adata.obs["lineage"]=="undefined")), "lineage_id"] = "Undefined"
adata.obs.loc[((adata.obs["cell_group"]=="2D")), "lineage_id"] = "2D EPSCs"

# Shows counts for all lineages.
adata.obs.lineage_id.value_counts()


adata.obs["name"] = "Undefined"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="Nat")), "name"] = "Epiblast"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="D5")), "name"] = "ELCs"
adata.obs.loc[((adata.obs["lineage"]=="epi") & (adata.obs["cell_group"]=="D6")), "name"] = "ELCs"

adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="Nat")), "name"] = "Hypoblast"
adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="D5")), "name"] = "HLCs"
adata.obs.loc[((adata.obs["lineage"]=="hypo") & (adata.obs["cell_group"]=="D6")), "name"] = "HLCs"

adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="Nat")), "name"] = "Trophectoderm"
adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="D5")), "name"] = "TLCs"
adata.obs.loc[((adata.obs["lineage"]=="te") & (adata.obs["cell_group"]=="D6")), "name"] = "TLCs"

adata.obs.loc[((adata.obs["lineage"]=="undefined")), "name"] = "Undefined"
adata.obs.loc[((adata.obs["cell_group"]=="2D")), "name"] = "2D EPSCs"


# Subset to look at each cell group
adata_2D = adata[adata.obs['cell_group']== '2D']
adata_pre = adata[adata.obs['cell_group']== 'D5']
adata_post = adata[adata.obs['cell_group']== 'D6']
adata_nat = adata[adata.obs['cell_group']== 'Nat']


adata_epsc_nat = adata[(adata.obs['cell_group'] != 'D5') &
                       (adata.obs['cell_group'] != 'D6')]

# Subset to look at each lineage
adata_lineage = adata[adata.obs['lineage']!= 'undefined']

# Subset excluding natural samples
adata_epscs = adata[adata.obs['cell_group'] != 'Nat']

# Subset with synthetic structures only
adata_synth = adata_epscs[adata_epscs.obs['cell_group'] != '2D']
adata_synth = adata_synth[adata_synth.obs['lineage_id'] != 'Undefined']
adata_synth = adata_synth[adata_synth.obs['lineage'] != 'undefined']

# Subset for EPI and ELCs
adata_epi = adata[adata.obs['lineage']== 'epi']
adata_epi = adata[(adata.obs['lineage_id']== 'Epiblast') |
                  (adata.obs['lineage_id']== 'D5 ELCs') |
                  (adata.obs['lineage_id']== 'D6 ELCs')]


# Subset for HYPO and HLCs
adata_hypo = adata[adata.obs['lineage']== 'hypo']
adata_hypo = adata[(adata.obs['lineage_id']== 'Hypoblast') |
                   (adata.obs['lineage_id']== 'D5 HLCs') |
                   (adata.obs['lineage_id']== 'D6 HLCs')]


# Subset for TE and TLCs
adata_te = adata[adata.obs['lineage']== 'te']
adata_te = adata[(adata.obs['lineage_id']== 'Trophectoderm') |
                 (adata.obs['lineage_id']== 'D5 TLCs') |
                 (adata.obs['lineage_id']== 'D6 TLCs')]


adata_name = adata[adata.obs['lineage'] != 'undefined']
# adata_name = adata_name[adata_name.obs['name'] != '2D EPSCs']

adata_lineage = adata[adata.obs['lineage']!= 'undefined']
adata_lineage = adata_lineage[adata_lineage.obs['name']!= '2D EPSCs']

palette_umap = {"Epiblast":'darkcyan',
                "Hypoblast":"darkgoldenrod",
                "Trophectoderm":"mediumvioletred",
                "Undefined":"whitesmoke",
                "D5 ELCs":'aquamarine',
                "D5 HLCs":"khaki",
                "D5 TLCs":"hotpink",
                "D6 ELCs":'turquoise',
                "D6 HLCs":"gold",
                "D6 TLCs":"deeppink",
                "nan":"whitesmoke",
                "epi":'lightseagreen',
                "hypo":"goldenrod",
                "te":"mediumvioletred",
                "undefined":"whitesmoke",
                "2D EPSCs":'royalblue',
                "2D":'royalblue',
                "ELCs":'lightseagreen',
                "HLCs":"goldenrod",
                "TLCs":"deeppink"}

# Assigns colors for UMAP.
colors_umap = {"Day 5 hEP-structures":'dodgerblue',
               "Day 6 hEP-structures":"hotpink",
               "Natural human embryo":"gold",
               "hEPSCs in 2D":"thistle"}


# Load in statistical analysis for HYPO related cells.
df_hypo_p = pd.read_excel("data/predict_data/Sozen_NC_2021/hEP-structures_MZG/data/Supp.Table4.xlsx", sheet_name=0, usecols=["gene name", "D5 Dunnet pval", "D6 Dunnet pval"])

# # Load in statistical analysis for EPI related cells.
df_epi_p = pd.read_excel("data/predict_data/Sozen_NC_2021/hEP-structures_MZG/data/Supp.Table4.xlsx", sheet_name=1, usecols=["gene name", "D5 Dunnet pval", "D6 Dunnet pval"])


# # Load in statistical analysis for TE related cells.
df_te_p = pd.read_excel("data/predict_data/Sozen_NC_2021/hEP-structures_MZG/data/Supp.Table4.xlsx", sheet_name=2, usecols=["gene name", "D5 Dunnet pval", "D6 Dunnet pval"])

# Function to add significance to violin plots.
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def add_significance(plot, x1, x2, values1, values2, p):
    y, h, col = max(max(values1),max(values2)) , 0.5, 'k'
    if y>30:
        h = 3
    if y<5:
        h=0.5
    if is_number(p):
        p=float(p)
        if p <= 0.0001 :
            TXT = "****"
        elif p <=  0.001:
            TXT = "***"
        elif p< 0.01:
            TXT = "**"
        elif p< 0.05:
            TXT = "*"
        else:
            TXT = "ns"
    else:
        TXT = "****"
    plot.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plot.text((x1+x2)*.5, y+h, TXT, ha='center', va='bottom', color=col)
    plot.set_ylim(top=y+h+1)
    if y>30:
        plot.set_ylim(top=y+h+3)
    if y<5:
        plot.set_ylim(top=y+h+0.5)


# Establish figure axis
fig5a, (ax0) = plt.subplots(1, 1, figsize=(10,8))

# Figure 5a for sample clustering UMAP
sc.pl.umap(adata, color="sample_group",
           add_outline=True,
           outline_width=[0.05, 0.01],
           size=100,
           legend_fontsize=20,
           palette = colors_umap,
           ax=ax0,
           wspace=5,
           title="  ",
           show=False)
ax0.set(xlabel=None, ylabel=None)

# Establish figure axis
fig5b, (ax1) = plt.subplots(1, 1, figsize=(10,8))

# Figure 5b for lineage scoring UMAP
sc.pl.umap(adata, color="name",
           add_outline=True,
           outline_width=[0.03, 0.01],
           size=75,
           legend_fontsize=20,
           cmap="bwr",
           palette = palette_umap,
           ax=ax1,
           wspace=5,
           title="  ",
           alpha=0.9,
           show=False)
ax1.set(xlabel=None, ylabel=None)

plt.rcParams.update({"font.size":50})
matplotlib.rcParams['axes.titlesize'] = 50

vmax_var = "p99"
vmin_var = "p10"

fig, axes= plt.subplots(2,3, figsize=(50,25), gridspec_kw={'wspace':0.1})
list_genes = ["GATA3", "TFAP2C", "KRT8",
              "ISL1", "TFAP2A", "KRT7"]

for i, marker in enumerate(list_genes):
    sc.pl.umap(adata,
               color=list_genes[i],
               add_outline=True,
               outline_width=[0.1, 0.01],
               vmax=vmax_var,
               vmin=vmin_var,
               ax=axes[int(i/3)][i%3],
               color_map="bwr",
               size=150,
               show=True,
               use_raw=False,
               title="  ")
    axes[int(i/3)][i%3].axis('off')
    axes[int(i/3)][i%3].collections[-1].colorbar.remove()

fig.text(0.13, 0.7, list_genes[0], fontname="sans-serif", fontstyle="italic")
fig.text(0.40, 0.7, list_genes[1], fontname="sans-serif", fontstyle="italic")
fig.text(0.67, 0.7, list_genes[2], fontname="sans-serif", fontstyle="italic")
fig.text(0.13, 0.27, list_genes[3], fontname="sans-serif", fontstyle="italic")
fig.text(0.40, 0.27, list_genes[4], fontname="sans-serif", fontstyle="italic")
fig.text(0.67, 0.27, list_genes[5], fontname="sans-serif", fontstyle="italic")

#adata.write("tmp_data/May2_2021/Sozen_2021_nc.h5ad")
pd.DataFrame(data=adata.raw.X.toarray(), index=adata.obs_names, columns=adata.raw.var_names).to_csv('data/predict_data/Sozen_NC_2021/rec.counts.t.csv')
adata.obs.to_csv("data/predict_data/Sozen_NC_2021/rec.meta.csv")
adata.obsm.to_df().to_csv("data/predict_data/Sozen_NC_2021/rec.cell.coord.csv")

#p2j rec.hEP-structures_MZG.py
#ipython3 nbconvert  rec.hEP-structures_MZG.ipynb --execute  --to html