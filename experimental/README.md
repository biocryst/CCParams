
### Anti-parallel dimers

Note that in this example all fragments are processed, producing an 'overall' parameterization, 
which is dominated by the rotation/translation and is roughly equivalent to the results on Fig. 2b,c, 
section 4 of the paper. To create a fine set of parameters, select fragments in a certain single 
register first. 

Also note that no changes are necessary in the code that transforms to and from the parameters. 
The fact that the two chains are anti-parallel is fully encoded in the rotation matrix.

1. Run *ProcessAntiDimerCoords.py*. Use the same helix template as for the parallel dimers.

Expected output:
```
CC reconstruction stats, RMSD:
min: 0.107829340322
max: 1.28285623225
mean: 0.262238338738
median: 0.230404333733
```


2. Run *AntiCCParamsStats.py*, it will produce file *pca_dimer_a_10.pkl*, 
which can be used for interactive visualisation.

Expected output:
```
Number of components: 10
Fraction of structures under 1A RMSD: 0.999452654625
CC reconstruction stats, RMSD:
min: 0.116282849824
max: 2.74863516976
mean: 0.318168269092
median: 0.287034712615
```

3. Launch interactive visualisation in PyMol, the same way as for the parallel dimer. 
Since the training set here contains all registers, the default mean position will have two helices overlapped.
The register here is controlled using first two parameters. 