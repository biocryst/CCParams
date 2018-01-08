# CCParams
Optimal data-driven parameterization of coiled coils

Requirements:  pyRMSD, scipy, numpy, sklearn, Biopython.

### cc_dataset_mmol_all
List of parallel dimeric CC structures exported from the CC+ database.

### ExtractDimerACoords.py
Extracting 15-residue CC fragments, here limited to the a-cluster.

### HelixTemplate.py
Parameterization of an a-helix.

### CCParamsLib.py
Parameterization of a CC fragment.

### ProcessDimerACoords.py
Generation of an extended parameterisation vector for every CC fragment.

### CCParamsStats.py
PCA of the extended dataset with varying number of parameters.

### CCParamsInteractive.py

Interactive visualisation of the CC parameters in PyMOL. Requires PyQT5, PyRosetta, and PyMOLPyRosettaServer link established in PyMOL (see PyRosetta documentation).

### helix_template.pkl and pca_dimer_a.pkl

Pickled PCA objects with parameterizations of an a-helical and a dimeric CC (a-cluster) fragments.
