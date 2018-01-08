# CCParams
Optimal data-driven parameterization of coiled coils

## Requirements
pyRMSD, scipy, numpy, sklearn, Biopython. 

## Description of files

### cc_dataset_mmol_all
List of parallel dimeric CC structures exported from the CC+ database, used in this study.

### ExtractDimerACoords.py
Extraction of 15-residue CC fragments, here limited to the a-cluster.

### HelixTemplate.py
Parameterization of an a-helix.

### CCParamsLib.py
Extended parameterization and reconstruction procedures.

### ProcessDimerACoords.py
Generation of an extended parameterisation vector for every CC fragment, preparing the dataset for PCA.

### CCParamsStats.py
PCA of the extended dataset and respective parameterizations with varying number of parameters.

### CCParamsInteractive.py

Interactive visualisation of the CC parameters in PyMOL. Requires PyQT5, PyRosetta, and the PyMOL-PyRosetta link must be established prior to running the app (see PyRosetta documentation).

### helix_template.pkl and pca_dimer_a.pkl

Pickled PCA objects with parameterizations of an a-helical and a dimeric CC (a-cluster) fragments.
