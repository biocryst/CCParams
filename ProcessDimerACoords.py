from CCParamsLib import *

cc_coords_all = cPickle.load(open('dimer_a_unique_0.2.pkl', "rb"))

n_cc = cc_coords_all.shape[0]
rmss = []
cc_params_all = []
for i in range(n_cc):
    cc = cc_coords_all[i]
    cc_params = cc2params(cc)
    cc_params_all.append(cc_params)
    cc_reconstructed = params2cc(cc_params)

    sup = SVDSuperimposer()
    sup.set(cc, cc_reconstructed)
    sup.run()
    rms = sup.get_rms()
    rmss.append(rms)

print "CC reconstruction stats, RMSD:\nmin: {}\nmax: {}\nmean: {}\nmedian: {}".format(np.min(rmss),np.max(rmss),np.mean(rmss),np.median(rmss))

cc_params_all = np.array(cc_params_all)
cPickle.dump((cc_params_all, cc_coords_all), open('dimer_a_params_coords.pkl', "wb"))
