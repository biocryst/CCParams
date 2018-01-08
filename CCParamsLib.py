import cPickle
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.linalg import sqrtm, inv


def sym(w):
    return w.dot(inv(sqrtm(w.T.dot(w))))

def params2cc(parameters, return_steps=False):
    steps = []

    h2_tran = parameters[-3:]

    h2_rot = np.reshape(parameters[-12:-3], (3,3))
    helical_params = np.reshape(parameters[:-12], (2,3))

    h2_rot = sym(h2_rot)

    h1_ref = du_mean_helix
    if return_steps:
        step_coords = np.append(h1_ref,h1_ref,axis=0)
        steps.append(step_coords)

    h1_dev, h2_dev = du_pca_helix.inverse_transform(helical_params)
    h1_dev = np.reshape(h1_dev, (h1_ref.shape[0],h1_ref.shape[1]))
    h2_dev = np.reshape(h2_dev, (h1_ref.shape[0],h1_ref.shape[1]))

    h1 = h1_ref+h1_dev
    h2 = h1_ref+h2_dev

    if return_steps:
        step_coords = np.append(h1,h2,axis=0)
        steps.append(step_coords)

    h2_new = np.dot(h2, h2_rot)
    if return_steps:
        step_coords = np.append(h1,h2_new,axis=0)
        steps.append(step_coords)

    h2_new = h2_new + h2_tran
    if return_steps:
        step_coords = np.append(h1,h2_new,axis=0)
        steps.append(step_coords)

    orig_coords = np.append(h1,h2_new,axis=0)

    if return_steps:
        return steps

    return orig_coords


def cc2params(coords):

    sup=SVDSuperimposer()

    n_atoms_mono = int(coords.shape[0]/2)

    h1 = coords[:n_atoms_mono]
    h2 = coords[n_atoms_mono:]
    h1_ref = du_mean_helix

    # align h1 and h2 with mean angles to the ref helix

    sup.set(h1_ref, h1)
    sup.run()
    h1_aligned_ref = sup.get_transformed()

    sup.set(h1_ref, h2)
    sup.run()
    h2_aligned_ref = sup.get_transformed()

    # estimate parameters from pca
    # center to h1_ref
    h1_aligned_ref = h1_aligned_ref - h1_ref
    h2_aligned_ref = h2_aligned_ref - h1_ref
    # unwrap
    h1_aligned_ref = np.reshape(h1_aligned_ref, (h1_aligned_ref.shape[0]*h1_aligned_ref.shape[1]))
    h2_aligned_ref = np.reshape(h2_aligned_ref, (h2_aligned_ref.shape[0]*h2_aligned_ref.shape[1]))
    # get params
    helical_params =du_pca_helix.transform([h1_aligned_ref, h2_aligned_ref])
    h1_aligned_ref, h2_aligned_ref = du_pca_helix.inverse_transform(helical_params)
    h1_aligned_ref = np.reshape(h1_aligned_ref, (h1_ref.shape[0],h1_ref.shape[1]))
    h2_aligned_ref = np.reshape(h2_aligned_ref, (h1_ref.shape[0],h1_ref.shape[1]))

    # construct ideal helices

    h1_transformed = h1_ref+h1_aligned_ref
    h2_transformed = h1_ref+h2_aligned_ref

    # adjust hi_helix1 and hi_helix2 by the parameters

    # align h1 to h1 ideal and transform all coords

    sup.set(h1, h1_transformed)
    sup.run()
    h1_ideal = sup.get_transformed()

    # align h2_ideal to h2

    sup.set(h2, h2_transformed)
    sup.run()
    h2_ideal = sup.get_transformed()

    coords_ideal = np.append(h1_ideal,h2_ideal,axis=0)

    # align ideal coords to the ref helix
    sup.set(h1_transformed,h1_ideal)
    sup.run()
    (rot_ref, tran_ref) = sup.get_rotran()

    coords_ideal = np.dot(coords_ideal,rot_ref) + tran_ref

    h1_new = coords_ideal[:n_atoms_mono]
    h2_new = coords_ideal[n_atoms_mono:]

    sup.set(h2_new, h2_transformed)
    sup.run()
    (rot2, tran2) = sup.get_rotran()
    h2_rot = rot2.flatten()
    helical_params = helical_params.flatten()
    transform_params = np.append(helical_params,np.append(h2_rot, tran2))

    return transform_params


du_pca_helix, du_mean_helix = cPickle.load(open('helix_template.pkl', "rb"))

