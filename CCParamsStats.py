from sklearn.decomposition import PCA
from scipy.optimize import minimize
from CCParamsLib import *


def min_rec(x):
    xa = np.abs(x)
    coefs_in = [xa[0],xa[1],np.abs(1-xa[0]-xa[1])]
    coefs_in = np.array(coefs_in)/np.sum(coefs_in)

    coef_hel = coefs_in[0]
    coef_rot = coefs_in[1]
    coef_angle = coefs_in[2]

    coefs = [coef_hel]*6+[coef_rot]*9+[coef_angle]*3

    training_set = cc_params_all*coefs
    n_components = 3
    pca_n = PCA(n_components=n_components)
    transformed = pca_n.fit_transform(training_set)
    reconstructed = pca_n.inverse_transform(transformed)
    reconstructed = reconstructed/coefs

    rmss = []
    sup = SVDSuperimposer()
    for i in test_inds:
        cc = cc_coords_all[i]

        cc_params = reconstructed[i]
        cc_reconstructed = params2cc(cc_params)

        sup.set(cc, cc_reconstructed)
        sup.run()
        rms = sup.get_rms()
        rmss.append(rms)
    print coefs_in, np.mean(rmss)
    return np.mean(rmss)

def optimise_params(transformed,coefs):

    def min_rec_cc(x):
        reconstructed = pca_n.inverse_transform([x])[0] / coefs
        cc_reconstructed = params2cc(reconstructed)

        sup = SVDSuperimposer()
        sup.set(orig_cc, cc_reconstructed)
        sup.run()
        rms = sup.get_rms()
        return rms

    rmss = []
    for ind_cc,short_params in enumerate(transformed):
        orig_cc = cc_coords_all[ind_cc]
        res = minimize(min_rec_cc,short_params)
        rms = res.fun
        rmss.append(rms)
    return rmss


cc_params_all, cc_coords_all = cPickle.load(open('dimer_a_params_coords.pkl', "rb"))

cc_coords_all = np.array(cc_coords_all)
cc_params_all = np.array(cc_params_all)

n_cc = cc_coords_all.shape[0]

# weight optimisation, uncomment to re-run
# test_inds = xrange(n_cc)
# x0 = [0.33,  0.33]
# res = minimize(min_rec, np.array(x0))
# x = res.x


x= [0.03438262,  0.79408506]
coefs_in = np.abs([x[0], x[1], 1 - x[0] - x[1]])
coefs_in = np.array(coefs_in) / np.sum(coefs_in)

coef_hel = coefs_in[0]
coef_rot = coefs_in[1]
coef_angle = coefs_in[2]

coefs = [coef_hel] * 6 + [coef_rot] * 9 + [coef_angle] * 3

training_set = cc_params_all * coefs

for n_components in range(1,10):
    pca_n = PCA(n_components=n_components)
    transformed = pca_n.fit_transform(training_set)
    cPickle.dump((pca_n,coefs), open('pca_dimer_a_'+str(n_components)+'.pkl', "wb"))

    rmss = optimise_params(transformed,coefs)

    rmss = np.array(rmss)

    cPickle.dump(rmss,open('cc_a_template_rmss'+str(n_components)+'.pkl', "wb"))
    print n_components
    print 'precent under 1',len(rmss[rmss<1])/len(rmss)
    print 'min,mean,median,max:',np.min(rmss),np.mean(rmss), np.median(rmss), np.max(rmss)