"""
04_scenicplus/04_pycistopic_lda.py

Builds a cisTopic object from the ArchR peak x cell matrix exported by
01_export_inputs.R and runs LDA topic modelling (30 topics) with
pycisTopic's collapsed Gibbs sampler (run_cgs_models). Note: Mallet is NOT
used, even though the output file name says "Mallet".

Methods: "the scATAC data was pre-processed with pycisTopic (n_topics = 30)"

Input : SCENIC_DIR/atac/pycisTopic_input/{peak_counts.mtx, peaks.tsv, cells.tsv, cell_metadata.csv}
Output: SCENIC_DIR/output/ATAC_06_24_25/ATAC_062425_Mallet_models_30topics.pkl

Next step: 05_prep_scenicplus_inputs.py (adds the model to the cisTopic object,
binarizes topics and attaches cell-type labels).
"""
import os
import pickle
import warnings

import pandas as pd
from scipy.io import mmread

import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models

warnings.simplefilter(action="ignore")
print("pycisTopic version:", pycisTopic.__version__)

# ---- Paths (override with environment variables) ------------------------------
projDir = os.environ.get("KRT17_SCENIC_DIR", "/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/scenic/")

inDir = os.path.join(projDir, "atac", "pycisTopic_input")
outDir = os.path.join(projDir, "output", "ATAC_06_24_25")
tmpDir = projDir
os.makedirs(outDir, exist_ok=True)

# ---- Create cisTopic object ----------------------------------------------------
print("Loading matrix and counts")
counts = mmread(os.path.join(inDir, "peak_counts.mtx")).tocsr()
peaks = pd.read_csv(os.path.join(inDir, "peaks.tsv"), header=None)[0].tolist()
cells = pd.read_csv(os.path.join(inDir, "cells.tsv"), header=None)[0].tolist()
cellmeta = pd.read_csv(os.path.join(inDir, "cell_metadata.csv"))

print("Making cistopic object")
cobj = create_cistopic_object(fragment_matrix=counts, cell_names=cells, region_names=peaks)
print(cobj)

# ---- LDA topic modelling (30 topics) ---------------------------------------------
print("Running topic modelling")
models = run_cgs_models(
    cobj,
    n_topics=[30],
    n_cpu=24,
    n_iter=300,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=tmpDir,
)

print("Saving models")
with open(os.path.join(outDir, "ATAC_062425_Mallet_models_30topics.pkl"), "wb") as f:
    pickle.dump(models, f)
