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

# Reproducing results of the paper

Here we focus only on the a-cluster parameterization. The procedure is, however, the same if one wishes 
to process fragments in other registers, or all at once. 

1. Extract PDB coordinates from *MMOL.tar.xz*. 
2. Run *ExtractDimerACoords.py*. This script will go through records in the *cc_dataset_mmol_all* file, 
locate all listed CC structures and extract regular fragments which start from the *a* heptad position. 
Expected script output:
```
Number of coiled coil records: 1866
Total number of fragments collected: 4118
Number of unique fragments: 2624
Number of unique fragments excluding outliers: 2538
```
Two files should appear, *dimer_a_all.pkl* and *dimer_a_unique_0.2.pkl*. 
The first one contains coordinates of all collected CC fragments from the a-cluster, 
the second one -- only non-redundant subset (in which no structure is similar with RMSD less than 0.2).

3. Run *HelixTemplate.py*. Note that in the paper we used coordinates of all fragments for this step, 
while here we take only the a-cluster as the input. This has no impact on the obtained parameterization and reconstruction
statistics (cf. Fig. 1d).

Expected script output:
```
Total number of helices: 8236
Number of unique helices: 2967
Helix reconstruction stats, RMSD:
min: 0.0779318211206
max: 0.792185424998
mean: 0.229877830959
median: 0.205617162334
```
