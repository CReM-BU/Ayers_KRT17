"""
04_scenicplus/05_prep_scenicplus_inputs.py

Prepares the two SCENIC+ inputs (L. Ayers):
  1. ATAC: rebuilds the cisTopic object, adds the 30-topic LDA model from
     04_pycistopic_lda.py, binarizes topics (otsu and top-3k region sets ->
     region_sets/, used by the Snakemake pipeline), and labels cells with
     their cell type (`ident`).
     Output: ATAC_062425_cistopicobj_30topics_idents.pkl
  2. RNA: sets AnnData .X to log-normalised counts.
     Output: rna/rna.fixed.h5ad

Requires: outputs of 01_export_inputs.R (pycisTopic_input/, rna/rna.h5ad)
          and 04_pycistopic_lda.py (ATAC_062425_Mallet_models_30topics.pkl)

Cell-type labels must be IDENTICAL in the ATAC (`ident` below) and RNA
(`ident` in rna.h5ad = subcluster_names) objects: SCENIC+ pairs ATAC and RNA
cells into metacells by this label (key_to_group_by: "ident" in
06_scenicplus_snakemake/config.yaml), and groups without a match are dropped.
"""
import logging
import sys
import os
from datetime import datetime

# Same directory as KRT17_SCENIC_DIR / SCENIC_DIR in config.R
projDir = os.environ.get('KRT17_SCENIC_DIR', '/Users/ljayers/Desktop/scenic').rstrip('/')

import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__

outDir = projDir + '/output/ATAC_06_24_25'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

tmpDir = projDir
# Create cisTopic object
from scipy.io import mmread
import pandas as pd
from pycisTopic.cistopic_class import *

counts = mmread(projDir+"/atac/pycisTopic_input/peak_counts.mtx")
peaks = pd.read_csv(projDir+"/atac/pycisTopic_input/peaks.tsv", header = None)
cells = pd.read_csv(projDir+"/atac/pycisTopic_input/cells.tsv", header = None)
counts = counts.tocsr()

# Make cistopic object
cells = cells[0].tolist()
peaks = peaks[0].tolist()
cobj = create_cistopic_object(fragment_matrix = counts, cell_names = cells, region_names = peaks)


# add cell metadata to cistopic object
cellmeta = pd.read_csv(projDir+"/atac/pycisTopic_input/cell_metadata.csv")
cellmeta['cell'] = cellmeta['cell'].astype(str) + '___cisTopic'
cellmeta = cellmeta.set_index('cell')
cobj.add_cell_data(cellmeta)

# load model from pycistopic, add to cistopic object
import pickle
infile = open(os.path.join(outDir, 'ATAC_062425_Mallet_models_30topics.pkl'), 'rb')  # written by 04_pycistopic_lda.py
model = pickle.load(infile)
infile.close()

model_select = evaluate_models(model,
                     select_model=30, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'loglikelihood'],
                     plot_metrics=False)
cobj.add_LDA_model(model_select)

with open(outDir + '/ATAC_cisTopicObject_30topics.pkl', 'wb') as f:
  pickle.dump(cobj, f)


## Topic binarization
if not os.path.exists(outDir+'/topic_binarization'):
    os.makedirs(outDir+'/topic_binarization')
from pycisTopic.topic_binarization import *
region_bin_topics_top_3k = binarize_topics(
    cobj, method='ntop', ntop = 3000,
    plot=False, num_columns=5)

region_bin_topics_otsu = binarize_topics(
    cobj, method='otsu',
    plot=False, num_columns=5)


binarized_cell_topic = binarize_topics(
    cobj,
    target='cell',
    method='li',
    plot=False,
    num_columns=5, nbins=100)

with open(outDir + '/test_binarized_cell_topic_30topics.pkl', 'wb') as f:
  pickle.dump(binarized_cell_topic, f)



os.makedirs(os.path.join(outDir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(outDir, "region_sets", "Topics_top_3k"), exist_ok = True)
#os.makedirs(os.path.join(outDir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.utils import region_names_to_coordinates
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(outDir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )



### Label cells with cell types (same mapping as 03_scATAC_analysis/02_integrate_cluster_annotate.R)
# Note: the heatmaps in results/ come from an earlier run that used older names
# (ABI1 = Pre KRT5-/KRT17+, ABI2 = KRT5-/KRT17+, TGFBUnk/TGFBUnk2, SecretoryLike,
# iBC1/iBC2) in both the ATAC and RNA objects.
file_path = outDir + '/ATAC_cisTopicObject_30topics.pkl'
import pickle
with open(file_path, 'rb') as file:
    cobj = pickle.load(file)

mapping = {
    'C1': 'iBC 2',
    'C2': 'iBC 1',
    'C3': 'PNEC-like',
    'C4': 'iAT1',
    'C5': 'SOX high TGFB-treated',   # was 'SOX high TGFB' -- must match the RNA label
    'C6': 'iAT2',
    'C7_1': 'KRT5-/KRT17+',
    'C7_2': 'Pre KRT5-/KRT17+',
    'C8': "Differentiating Airway"
}

cobj.cell_data['ident'] = cobj.cell_data['subclusters'].replace(mapping)


# Analyze topics
from pycisTopic.topic_qc import topic_annotation
topic_annot = topic_annotation(cobj, annot_var='ident', binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)
topic_annot


### Save annotated cistopic object
with open(outDir + '/ATAC_062425_cistopicobj_30topics_idents.pkl', 'wb') as f:
  pickle.dump(cobj, f)



### Prep RNA .h5ad file (Set adata.X to the logcounts layer)
import anndata as ad
import numpy as np

adata = ad.read_h5ad(os.path.join(projDir, "rna", "rna.h5ad"))

# save original X matrix
adata.raw = adata.copy()

# set main expression matrix to logcounts
adata.X = adata.layers["logcounts"].copy()

# save fixed object
adata.write(
    os.path.join(projDir, "rna", "rna.fixed.h5ad"),
    compression="gzip"
)



