import cPickle
import numpy as np
from pyRMSD import RMSDCalculator
from Bio.SVDSuperimposer import SVDSuperimposer
from sklearn.decomposition import PCA


cc_coords_all = cPickle.load(open('dimer_a_all.pkl', "rb"))

hel_coords_all = np.reshape(cc_coords_all,(cc_coords_all.shape[0]*2,int(cc_coords_all.shape[1]/2),cc_coords_all.shape[2]))
n_helices = hel_coords_all.shape[0]

print "Total number of helices:", n_helices

# filter out structurally redundant a-helices
threshold = 0.1
h_unique = [hel_coords_all[-1]]
global_ind = 0
unique_ind = 0
for ah in reversed(hel_coords_all[:-1]):
    h_test = [ah] + h_unique
    h_test = np.array(h_test)
    calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", h_test)
    dist = calculator.oneVsFollowing(0)
    global_ind += 1
    if np.min(dist) > threshold:
        h_unique.append(ah)
        unique_ind += 1

h_unique = np.array(h_unique)
n_helices = h_unique.shape[0]
print "Number of unique helices:", n_helices

# superpose all
calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", h_unique)
calculator.iterativeSuperposition()

# center data for PCA
mean_helix = np.mean(h_unique,0)
straight_helices_al = h_unique - mean_helix
mean_helix -= np.mean(mean_helix,0)

XYZ_flat = np.reshape(straight_helices_al, (n_helices, straight_helices_al.shape[1]*straight_helices_al.shape[2]))

n_components = 3
pca_n = PCA(n_components=n_components)
transformed = pca_n.fit_transform(XYZ_flat)

reconstructed_flat = pca_n.inverse_transform(transformed)
reconstructed = np.reshape(reconstructed_flat, straight_helices_al.shape)+mean_helix

# mean helix and 3 PCs make a helix template
cPickle.dump((pca_n,mean_helix), open('helix_template.pkl', "wb"))

# check that we can restore the geometry of all input helices with reasonable RMSD
rmss = []
for i in range(n_helices):
    sup=SVDSuperimposer()
    sup.set(straight_helices_al[i]+mean_helix,reconstructed[i])
    sup.run()
    rms = sup.get_rms()
    rmss.append(rms)

print "Helix reconstruction stats, RMSD:\nmin: {}\nmax: {}\nmean: {}\nmedian: {}".format(np.min(rmss),np.max(rmss),np.mean(rmss),np.median(rmss))
