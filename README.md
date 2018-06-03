# CCParams
Optimal data-driven parameterization of coiled coils

## Requirements
pyRMSD, scipy, numpy, sklearn, Biopython, PyRosetta (optional). 

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
statistics.

Expected script output (cf. Fig. 1d):
```
Total number of helices: 8236
Number of unique helices: 2967
Helix reconstruction stats, RMSD:
min: 0.0779318211206
max: 0.792185424998
mean: 0.229877830959
median: 0.205617162334
```

4. Run *ProcessDimerACoords.py*. This step will convert CC fragment coordinates from file 
*dimer_a_unique_0.2.pkl* into 'extended' parameterization and test the reconstruction accuracy.
The results will be stored in the file *dimer_a_params_coords.pkl*.

Expected output (cf. Fig 3c):
```
CC reconstruction stats, RMSD:
min: 0.0899008249849
max: 0.680252038429
mean: 0.232043925963
median: 0.212188994827
```

5. Run *CCParamsStats.py*. This script will go through CC parameterizations with the number of parameters 
varying from 1 to 10 and output reconstruction statistics. Every reconstruction's parameters are optimised to minimize RMSD
to the target structures, so it'll take a while. 

Expected output, the first and the last batches are shown (cf. Fig 3c):
```
Number of components: 1
Fraction of structures under 1A RMSD: 0.982269503546
CC reconstruction stats, RMSD:
min: 0.190088304232
max: 1.10620154112
mean: 0.517878054167
median: 0.48013578893
...

Number of components: 10
Fraction of structures under 1A RMSD: 1.0
CC reconstruction stats, RMSD:
min: 0.0968117019343
max: 0.685145833173
mean: 0.250808975295
median: 0.231687949811
```

6. Launch PyMOL, Open File-> Run Script..., find and select *PyMOL-RosettaServer.py* ([PyRosetta](http://www.pyrosetta.org/dow) is required) 
Run *CCParamsInteractive.py*. If everything went well, you should be able to 
modify parameters of CCs (first 10 are shown) and see the 
geometry changing in PyMOL in real time. cf. Fig. 3d-f