import logging
import sys
import os
from datetime import datetime


projDir = '/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/scenic/'


import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__

os.environ['MALLET_MEMORY'] = '180G'

outDir = projDir + 'output/ATAC_06_24_25'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

tmpDir = projDir
print('Loading matrix and counts')
# Create cisTopic object
from scipy.io import mmread
import pandas as pd
from pycisTopic.cistopic_class import *


counts = mmread(projDir+"atac/pycisTopic_input/peak_counts.mtx")
peaks = pd.read_csv(projDir+"atac/pycisTopic_input/peaks.tsv", header = None)
cells = pd.read_csv(projDir+"atac/pycisTopic_input/cells.tsv", header = None)
cellmeta = pd.read_csv(projDir+"atac/pycisTopic_input/cell_metadata.csv")
counts = counts.tocsr()

cells = cells[0].tolist()
peaks = peaks[0].tolist()
print('Making cistopic object')
cobj = create_cistopic_object(fragment_matrix = counts, cell_names = cells, region_names = peaks)

print(cobj)

import os

print('running topic modeling')
path_to_mallet_binary="Mallet-202108/bin/mallet"
# Run models
from pycisTopic.lda_models import run_cgs_models
models=run_cgs_models(
                    cobj,
                    n_topics=[30],
                    n_cpu=24,
                    n_iter=300, 
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=tmpDir)

print('Saving models')

# Save
import pickle 
with open(outDir+'/ATAC_062425_Mallet_models_30topics.pkl', 'wb') as f:
  pickle.dump(models, f)