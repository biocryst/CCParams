import cPickle
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.linalg import sqrtm, inv


def sym(w):
    return w.dot(inv(sqrtm(w.T.dot(w))))

def params2cc(parameters):

    h2_tran = parameters[-15:-12]
    h2_rot = np.reshape(parameters[-24:-15], (3,3))

    h3_tran = parameters[-3:]
    h3_rot = np.reshape(parameters[-12:-3], (3,3))

    helical_params = np.reshape(parameters[:-24], (3,3))

    h2_rot = sym(h2_rot)
    h3_rot = sym(h3_rot)

    h1_ref = du_mean_helix

    h1_dev, h2_dev, h3_dev = du_pca_helix.inverse_transform(helical_params)
    h1_dev = np.reshape(h1_dev, (h1_ref.shape[0],h1_ref.shape[1]))
    h2_dev = np.reshape(h2_dev, (h1_ref.shape[0],h1_ref.shape[1]))
    h3_dev = np.reshape(h3_dev, (h1_ref.shape[0],h1_ref.shape[1]))

    h1 = h1_ref+h1_dev
    h2 = h1_ref+h2_dev
    h3 = h1_ref+h3_dev


    h2_new = np.dot(h2, h2_rot)
    h2_new = h2_new + h2_tran
    h3_new = np.dot(h3, h3_rot)
    h3_new = h3_new + h3_tran

    orig_coords = np.append(h1,np.append(h2_new,h3_new,axis=0),axis=0)

    return orig_coords


def cc2params(coords):

    sup=SVDSuperimposer()

    n_atoms_mono = int(coords.shape[0]/3)

    h1 = coords[:n_atoms_mono]
    h2 = coords[n_atoms_mono:2*n_atoms_mono]
    h3 = coords[2*n_atoms_mono:]

    h1_ref = du_mean_helix

    # align h1 and h2 with mean angles to the ref helix

    sup.set(h1_ref, h1)
    sup.run()
    h1_aligned_ref = sup.get_transformed()

    sup.set(h1_ref, h2)
    sup.run()
    h2_aligned_ref = sup.get_transformed()

    sup.set(h1_ref, h3)
    sup.run()
    h3_aligned_ref = sup.get_transformed()

    # estimate parameters from pca
    # center to h1_ref
    h1_aligned_ref = h1_aligned_ref - h1_ref
    h2_aligned_ref = h2_aligned_ref - h1_ref
    h3_aligned_ref = h3_aligned_ref - h1_ref
    # unwrap
    h1_aligned_ref = np.reshape(h1_aligned_ref, (h1_aligned_ref.shape[0]*h1_aligned_ref.shape[1]))
    h2_aligned_ref = np.reshape(h2_aligned_ref, (h2_aligned_ref.shape[0]*h2_aligned_ref.shape[1]))
    h3_aligned_ref = np.reshape(h3_aligned_ref, (h3_aligned_ref.shape[0]*h3_aligned_ref.shape[1]))

    # get params
    helical_params =du_pca_helix.transform([h1_aligned_ref, h2_aligned_ref,h3_aligned_ref])
    h1_aligned_ref, h2_aligned_ref, h3_aligned_ref = du_pca_helix.inverse_transform(helical_params)
    h1_aligned_ref = np.reshape(h1_aligned_ref, (h1_ref.shape[0],h1_ref.shape[1]))
    h2_aligned_ref = np.reshape(h2_aligned_ref, (h1_ref.shape[0],h1_ref.shape[1]))
    h3_aligned_ref = np.reshape(h3_aligned_ref, (h1_ref.shape[0],h1_ref.shape[1]))

    # construct ideal helices

    h1_transformed = h1_ref+h1_aligned_ref
    h2_transformed = h1_ref+h2_aligned_ref
    h3_transformed = h1_ref+h3_aligned_ref

    # adjust hi_helix1 and hi_helix2 by the parameters

    # align h1 to h1 ideal and transform all coords

    sup.set(h1, h1_transformed)
    sup.run()
    h1_ideal = sup.get_transformed()

    # align h2_ideal to h2

    sup.set(h2, h2_transformed)
    sup.run()
    h2_ideal = sup.get_transformed()

    sup.set(h3, h3_transformed)
    sup.run()
    h3_ideal = sup.get_transformed()

    coords_ideal = np.append(h1_ideal,h2_ideal,axis=0)
    coords_ideal = np.append(coords_ideal,h3_ideal,axis=0)

    # align ideal coords to the ref helix
    sup.set(h1_transformed,h1_ideal)
    sup.run()
    (rot_ref, tran_ref) = sup.get_rotran()

    coords_ideal = np.dot(coords_ideal,rot_ref) + tran_ref

    h1_new = coords_ideal[:n_atoms_mono]
    h2_new = coords_ideal[n_atoms_mono:2*n_atoms_mono]
    h3_new = coords_ideal[2*n_atoms_mono:]

    sup.set(h2_new, h2_transformed)
    sup.run()
    (rot2, tran2) = sup.get_rotran()
    h2_rot = rot2.flatten()

    sup.set(h3_new, h3_transformed)
    sup.run()
    (rot3, tran3) = sup.get_rotran()
    h3_rot = rot3.flatten()

    helical_params = helical_params.flatten()
    h2_position_params = np.append(h2_rot, tran2)
    h3_position_params = np.append(h3_rot, tran3)
    transform_params = np.append(helical_params,np.append(h2_position_params,h3_position_params))

    return transform_params


du_pca_helix, du_mean_helix = cPickle.load(open('helix_template.pkl', "rb"))

