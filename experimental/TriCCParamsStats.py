from TriCCParamsLib import *
from sklearn.decomposition import PCA

cc_params_all, cc_coords_all= cPickle.load(open('trimer_params_coords.pkl', "rb"))

x0 = [0.01775954,  0.83417166,  0.1480688]
coef_hel = x0[0]
coef_rot = x0[1]
coef_tran = x0[2]
coefs = [coef_hel] * 9 + [coef_rot] * 9 + [coef_tran] * 3 + [coef_rot] * 9 + [coef_tran] * 3

training_set = cc_params_all * coefs

n_components = 10
pca_n = PCA(n_components=n_components)
transformed = pca_n.fit_transform(training_set)
cPickle.dump((pca_n,coefs), open('pca_trimer_'+str(n_components)+'.pkl', "wb"))

reconstructed = pca_n.inverse_transform(transformed) / coefs

n_cc = cc_coords_all.shape[0]
rmss = []
for i in range(n_cc):
    cc = cc_coords_all[i]

    cc_params = reconstructed[i]
    cc_reconstructed = params2cc(cc_params)

    sup = SVDSuperimposer()
    sup.set(cc, cc_reconstructed)
    sup.run()
    rms = sup.get_rms()
    rmss.append(rms)

rmss = np.array(rmss)
cPickle.dump(rmss,open('cc_template_tri_rmss.pkl', "wb"))

print "\nNumber of components:",n_components
print 'Fraction of structures under 1A RMSD:',float(len(rmss[rmss<1]))/float(len(rmss))
print "CC reconstruction stats, RMSD:\nmin: {}\nmax: {}\nmean: {}\nmedian: {}".format(np.min(rmss),np.max(rmss),np.mean(rmss),np.median(rmss))
