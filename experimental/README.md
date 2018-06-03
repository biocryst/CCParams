## Applications to other CC oligomeric states

Disclaimer: the code/data here are provided for demonstration purposes only. 
The data has not been cleaned, the code has neither been optimized nor validated 
on a sufficiently large set of experimental structures. 

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


2. Run *AntiCCParamsStats.py*. It will produce the file *pca_dimer_a_10.pkl*, 
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
Since the training set here contains all registers, it is controlled with the first two parameters.
The default mean position will have the two helices overlapped and PyMOL renderer somewhat confused. 
As a workaround, first adjust the sliders to distance the helices, delete the object in PyMOL and 
move the sliders again. 


### Trimers

Minor changes are required for the parametrisation and visualisation code, since we are working with three, 
rather than two helices. The adjusted files are provided in full as *TriCCParamsLib.py* and *TriCCParamsInteractive.py*.

1. Run *ProcessTrimerCoords.py*. Use the same helix template as for the parallel dimers.

Expected output:
```
CC reconstruction stats, RMSD:
min: 0.0941925679087
max: 0.648050494848
mean: 0.210258604266
median: 0.194409653232
```

2. Run *TriCCParamsStats.py*, it will produce file *pca_trimer_10.pkl*, 
   which can be used for interactive visualisation.

Expected output:
```
Number of components: 10
Fraction of structures under 1A RMSD: 0.999903138318
CC reconstruction stats, RMSD:
min: 0.127235329958
max: 1.06181204303
mean: 0.326316370077
median: 0.300639477533
```

3. Launch interactive visualisation with PyMol and *TriCCParamsInteractive.py*. 
Here again the first two parameters will control register and the default position 
will have three helices overlapped.