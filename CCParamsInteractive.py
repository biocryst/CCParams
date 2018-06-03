from __future__ import division
from PyQt5 import QtCore, QtWidgets
import sys, time
import pyrosetta
import pyrosetta.rosetta as rosetta
from CCParamsLib import *
from Bio.SVDSuperimposer import SVDSuperimposer


class AppForm(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setWindowTitle('CC parametrization')

        self.param_levels = []
        self.boxes = []

        self.pca_kr,self.coefs = cPickle.load(open('pca_dimer_a_10.pkl', "rb"))
        self.n_components = self.pca_kr.n_components
        self.sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        self.stepwise = False

        pyrosetta.init()
        self.pose = pyrosetta.Pose()
        pyrosetta.make_pose_from_sequence(self.pose, self.sequence, 'fa_standard',auto_termini=False)

        coord_params = np.zeros((1,self.n_components))
        coord_params = self.pca_kr.inverse_transform(coord_params)[0]/self.coefs
        self.mean_coords = params2cc(coord_params)
        self.mean_coords -= self.mean_coords.mean()
        self.sup = SVDSuperimposer()

        dummy = rosetta.numeric.xyzVector_double_t()
        dummy0 = rosetta.numeric.xyzVector_double_t()
        dummy0.x = 0
        dummy0.y = 0
        dummy0.z = 0
        ind = -1
        for r in range(self.pose.total_residue()):
            res = self.pose.residue(r+1)
            for a in range(res.natoms()):
                if self.pose.residue(r+1).atom_is_hydrogen(a+1) or self.pose.residue(r+1).atom_name(a+1) == ' OXT':
                    self.pose.residue(r + 1).set_xyz(a + 1, dummy0)
                    continue
                ind += 1
                v = self.mean_coords[ind]
                dummy.x = v[0]
                dummy.y = v[1]
                dummy.z = v[2]
                self.pose.residue(r+1).set_xyz(a+1,dummy)

        self.pmm = pyrosetta.PyMOLMover()
        self.pmm.keep_history(False)
        self.pmm.apply(self.pose)

        self.level_mult = [100]*self.n_components
        for i_par in range(self.n_components):

            box = QtWidgets.QHBoxLayout()
            label = QtWidgets.QLabel('X'+str(i_par+1))
            level = QtWidgets.QSlider(QtCore.Qt.Horizontal)

            level.setRange(-200, 200)

            level.setValue(0)
            level.setTracking(True)
            level.setEnabled(True)
            level.setTickPosition(QtWidgets.QSlider.TicksBelow)

            level.valueChanged.connect(self.on_change_level)

            box.addWidget(label)
            box.addWidget(level)
            self.boxes.append(box)
            self.param_levels.append(level)

        bt_reset = QtWidgets.QPushButton("Reset all")
        bt_reset.clicked.connect(self.reset_levels)

        buttons_box = QtWidgets.QHBoxLayout()
        buttons_box.addWidget(bt_reset)
        buttons_box.addStretch(1)

        top_box = QtWidgets.QVBoxLayout()

        for box in self.boxes:
            top_box.addLayout(box)

        top_box.addStretch(1)
        top_box.addLayout(buttons_box)

        parentWidget = QtWidgets.QWidget()
        parentWidget.setLayout(top_box)
        self.setCentralWidget(parentWidget)
        self.resize(400,600)

    def center_coords(self,coords):
        ref_cas = self.mean_coords[1::5]
        al_cas = coords[1::5]

        self.sup.set(ref_cas, al_cas)
        self.sup.run()
        rot,tran = self.sup.get_rotran()
        return np.dot(coords,rot)+tran

    def reset_levels(self):
        for level in self.param_levels:
            level.setValue(0)

    def on_change_level(self):
        coord_params = np.zeros((1,self.n_components))
        for p_ind,level in enumerate(self.param_levels):
            coord_params[0,p_ind] = level.value()/self.level_mult[p_ind]

        coord_params = self.pca_kr.inverse_transform(coord_params)[0]/self.coefs

        if self.stepwise:
            coords_list = params2cc(coord_params,True)
        else:
            coords_list = [params2cc(coord_params)]

        for coords in coords_list:
            coords = self.center_coords(coords)
            dummy = rosetta.numeric.xyzVector_double_t()
            ind = -1
            for r in range(self.pose.total_residue()):
                res = self.pose.residue(r+1)
                for a in range(res.natoms()):
                    if self.pose.residue(r+1).atom_is_hydrogen(a+1) or self.pose.residue(r+1).atom_name(a+1) == ' OXT':
                        continue
                    ind += 1
                    v = coords[ind]
                    dummy.x = v[0]
                    dummy.y = v[1]
                    dummy.z = v[2]
                    self.pose.residue(r+1).set_xyz(a+1,dummy)

            self.pmm.apply(self.pose)
            if self.stepwise:
                time.sleep(0.1)

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()
