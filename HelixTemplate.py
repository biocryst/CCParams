import random
import cPickle
import numpy as np
from pyRMSD import RMSDCalculator
from Bio.SVDSuperimposer import SVDSuperimposer
from sklearn.decomposition import PCA


cc_coords_all= cPickle.load(open('dimer_all.pkl', "rb"))
cc_coords_all = np.array(cc_coords_all)

hel_coords_all = np.reshape(cc_coords_all,(cc_coords_all.shape[0]*2,int(cc_coords_all.shape[1]/2),cc_coords_all.shape[2]))
n_helices = hel_coords_all.shape[0]

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
        if random.randint(1, 100) == 5:
            print global_ind, unique_ind

print "final:", global_ind, unique_ind

pca_n, mean_helix = cPickle.load(open('helix_template.pkl', "rb"))

calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", np.append([mean_helix], hel_coords_all,axis=0))
dist,straight_helices_al = calculator.oneVsTheOthers(0,True)

straight_helices_al = straight_helices_al[1:]
n_helices = straight_helices_al.shape[0]

calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", straight_helices_al)
calculator.iterativeSuperposition()

mean_helix2 = np.mean(straight_helices_al,0)

straight_helices_al = straight_helices_al - mean_helix2
mean_helix_center = np.mean(mean_helix2,0)
mean_helix2 = mean_helix2 - mean_helix_center

XYZ_flat = np.reshape(straight_helices_al, (n_helices, straight_helices_al.shape[1]*straight_helices_al.shape[2]))

n_components = 3
pca_n = PCA(n_components=n_components)
transformed = pca_n.fit_transform(XYZ_flat)

reconstructed_flat = pca_n.inverse_transform(transformed)
reconstructed = np.reshape(reconstructed_flat, straight_helices_al.shape)+mean_helix2

cPickle.dump((pca_n,mean_helix2), open('helix_template.pkl', "wb"))

rmss = []
for i in range(n_helices):
    sup=SVDSuperimposer()
    sup.set(straight_helices_al[i]+mean_helix2,reconstructed[i])
    sup.run()
    rms = sup.get_rms()
    rmss.append(rms)

print n_components, np.mean(rmss), np.median(rmss), np.max(rmss)
