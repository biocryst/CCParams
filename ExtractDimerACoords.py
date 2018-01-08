import os
from Bio.SeqUtils import seq1
import random
import warnings
import numpy as np
import cPickle
from Bio.PDB.Vector import rotaxis
from pyRMSD import RMSDCalculator

def get_gly_cb_coords(residue):
    try:
        n_v=residue["N"].get_vector()
        c_v=residue["C"].get_vector()
        ca_v=residue["CA"].get_vector()
    except:
        return None
    n_v=n_v-ca_v
    c_v=c_v-ca_v
    rot=rotaxis(-np.pi*120.0/180.0, c_v)
    cb_at_origin_v=n_v.left_multiply(rot)
    cb_v=cb_at_origin_v+ca_v
    return cb_v.get_array()

def split_cc(model, selection, n_res_split):
    (pdb_id, chains, start_res, end_res, sequences, registers) = selection
    chain1 = model[chains[0]]
    chain2 = model[chains[1]]

    n_res = end_res[0] - start_res[0] + 1
    backbone_ids = ['N', 'CA', 'C', 'O', 'CB']
    n_atoms_mono = n_res * len(backbone_ids)
    coords_all = np.zeros((n_atoms_mono * 2, 3))
    coord_ind = 0
    for res1_id, res2_id, aa1, aa2 in zip(range(start_res[0], end_res[0] + 1), range(start_res[1], end_res[1] + 1), sequences[0], sequences[1]):
        try:
            res1 = chain1[res1_id]
        except:
            res1 = chain1[('H_MSE', res1_id, ' ')]
        try:
            res2 = chain2[res2_id]
        except:
            res2 = chain2[('H_MSE', res2_id, ' ')]
        if aa1 != 'X':
            assert aa1 == seq1(res1.resname)
        if aa2 != 'X':
            assert aa2 == seq1(res2.resname)
        for backbone_id in backbone_ids:
            try:
                coords_all[coord_ind, :] = res1[backbone_id].coord
            except:
                assert backbone_id == 'CB'
                coords_all[coord_ind, :] = get_gly_cb_coords(res1)
            try:
                coords_all[n_atoms_mono + coord_ind, :] = res2[backbone_id].coord
            except:
                assert backbone_id == 'CB'
                coords_all[n_atoms_mono + coord_ind, :] = get_gly_cb_coords(res2)
            coord_ind += 1

    coiledcoils = []
    n_atoms_mono_split = n_res_split * len(backbone_ids)
    ref_reg = 'abcdefgabcdefga'
    for start_ind in range(0, n_res - n_res_split + 1):
        reg1 = registers[0][start_ind:start_ind+n_res_split]
        reg2 = registers[1][start_ind:start_ind+n_res_split]
        if reg1 != ref_reg or reg2 != ref_reg:
            continue
        start_ind_atom = start_ind * len(backbone_ids)
        coords1 = coords_all[start_ind_atom:start_ind_atom + n_atoms_mono_split, :]
        coords2 = coords_all[n_atoms_mono + start_ind_atom:n_atoms_mono + start_ind_atom + n_atoms_mono_split, :]
        mean1 = np.mean(coords1, 0)
        mean2 = np.mean(coords2, 0)
        dist = np.linalg.norm(mean1 - mean2)
        ca_dists = np.linalg.norm(coords1[1::5]-coords2[1::5],axis=1)
        if dist > 18 or np.max(ca_dists) > 18:
            print "wrong CC segment in", selection
            continue

        coords = np.append(coords1, coords2, axis=0)
        coiledcoils.append(coords)

    return coiledcoils


f = open('cc_dataset_mmol_all').read().splitlines()

records = [line.split() for line in f]
ccoils = []

for c1,c2 in zip(records[0::2], records[1::2]):
    assert c1[0] == c2[0]
    filename = c1[0]
    chains = (c1[1],c2[1])
    if len(c1[4]) != len(c2[4]):
        continue
    sequence1 = c1[4]
    sequence2 = c2[4]
    register1 = c1[5]
    register2 = c2[5]
    start_res = (int(c1[2]),int(c2[2]))
    end_res = (int(c1[3]),int(c2[3]))
    sequence = (sequence1,sequence2)
    register = (register1,register2)

    ccoils.append((filename, chains, start_res, end_res, sequence,register))
    ccoils.append((filename, chains[::-1], start_res[::-1], end_res[::-1], sequence[::-1],register[::-1]))

print len(ccoils)


pdb_path = "../CCFold/MMOL"
warnings.filterwarnings("ignore")

n_window = 15
cc_coords_all = []

for item in ccoils:
    (filename, chains, start_res, end_res, sequence,register) = item
    full_filename = os.path.join(pdb_path,filename)
    if not os.path.exists(full_filename):
        print 'missing ', full_filename
        continue
    structure = PDB.PDBParser().get_structure(filename, full_filename)
    model = structure[0]
    try:
        coiled_coils, cc_ca_dists = split_cc(model, item, n_window)
    except:
        print 'problem with', item
        print 'actual chains: ', structure.child_list, model.child_list
        continue
    for cc in coiled_coils:
        cc_coords_all.append(cc)

cPickle.dump(cc_coords_all, open('dimer_a_all.pkl', "wb"))

threshold = 0.2
cc_unique = [cc_coords_all[-1]]
global_ind = 0
unique_ind = 0
for cc in reversed(cc_coords_all[:-1]):
    cc_test = [cc]+cc_unique
    cc_test = np.array(cc_test)
    calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", cc_test)
    dist = calculator.oneVsFollowing(0)
    global_ind+=1
    if np.min(dist) > threshold:
        cc_unique.append(cc)
        unique_ind+=1

print "final:",global_ind,unique_ind

cc_unique = np.array(cc_unique)
calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", cc_unique)
calculator.iterativeSuperposition()
mean_cc = np.mean(cc_unique,0)

calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", np.append([mean_cc], cc_unique,axis=0))
dist = calculator.oneVsFollowing(0)

med_dist = np.median(dist)
mad = np.median(np.abs(dist-med_dist))
std_est = 1.4856*mad

cutoff = np.mean(dist)+3*np.std(dist)
idx = dist < cutoff
cc_filt = cc_unique[idx]
cPickle.dump(cc_filt, open('dimer_a_unique_0.2.pkl', "wb"))
print std_est
