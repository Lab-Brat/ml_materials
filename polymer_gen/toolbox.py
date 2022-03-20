import numpy as np

class tools():
    def __init__(self):
        self.pi = 3.1415926

        # polymer chain specific settings
        self.cca = ((180-109.5)*self.pi/180)/2
        self.cha = self.pi/3
        self.R, self.RT, self.RH, self.chainvec = self.get_matricies(self.cca, self.cha)
        self.bond_length = {'ccl': 1.33, 'chl': 0.99,
                            'nnl': 1.24, 'nhl': 1.01,
                            'cnl': 1.28, 'ncl': 1.28}

    def get_matricies(self, cca, cha):
        R = np.array([[np.cos(cca), np.sin(cca), 0],
                      [np.sin(-cca), np.cos(cca), 0],
                      [0, 0, 1]]).transpose()

        RT = np.array([[0, 1, 0],
                       [-1, 0, 0],
                       [0, 0, 1]]).transpose()

        RH = np.array([[1, 0, 0],
                       [0, np.cos(cha), np.sin(cha)],
                       [0, np.sin(-cha), np.cos(cha)]])

        chainvec = np.array([1, 0, 0])

        return R, RT, RH, chainvec

    def transpose(self):
        self.R = self.R.transpose()
        self.RT = self.RT.transpose()

    def gen_atom_init(self, bonds, cx):
        l1 = f"{bonds[0].lower()}l"
        l2 = f"{bonds[0].lower()}l"

        A1 = np.array([cx[0],
                       cx[1]+self.bond_length[l1]*np.cos(self.cha),
                       cx[2]+self.bond_length[l1]*np.sin(self.cha)])
        A2 = np.array([cx[0],
                       cx[1]+self.bond_length[l2]*np.cos(self.cha),
                       cx[2]+self.bond_length[l2]*np.sin(-self.cha)])

        return A1, A2

    def gen_atom_main(self, bond_main, prev_atom):
        l1 = f"{bond_main.lower()}l"
        AM = prev_atom + (self.bond_length[l1] * self.R.dot(self.chainvec))

        return AM

    def gen_atom_chain(self, bonds, last_a):
        l1 = f"{bonds[0].lower()}l"
        l2 = f"{bonds[0].lower()}l"

        A1 = last_a + self.bond_length[l1] * \
            (self.RH.dot(self.RT).dot(self.chainvec))
        A2 = last_a + self.bond_length[l2] * \
            (self.RH.transpose().dot(self.RT).dot(self.chainvec))

        return A1, A2

    def type_atom_size(self, t):
        if t == 'C':
            return 12.011
        elif t == 'N':
            return 14.006
        elif t == 'H':
            return 1.008
        else:
            print("[NO SUCH ELEMENT IN THE DATABASE]")

    def type_atom_ind(self, t):
        atom_ind = [1, 2, 3]

        if t == 'C':
            return atom_ind[0]
        elif t == 'N':
            return atom_ind[1]
        elif t == 'H':
            return atom_ind[2]
        else:
            print("[NO SUCH ELEMENT IN THE DATABASE]")

    def write_data(self, atom_types, atoms, bonds, angles, dihedrals, box_dim):
        with open('/home/labbrat/gitlab/polymer_gen/carbon_chain.data', 'w') as f:
            f.write(f"# Model for PE for {len(atoms)} atoms\n\n")
            f.write(f"     {len(atoms)}     atoms\n")
            f.write(f"     {len(bonds)}     bonds\n")
            f.write(f"     {len(angles)}     angles\n")
            f.write(f"     {len(dihedrals)}     dihedrals\n\n")
            f.write(f"     {len(set([i[1] for i in atoms]))}     atom types\n")
            f.write(f"     {len(set([i[1] for i in bonds]))}     bond types\n")
            f.write(f"     {len(set([i[1] for i in angles]))}     angle types\n")
            f.write(f"     {len(set([i[1] for i in dihedrals]))}     dihedral types\n\n")
            f.write(f"     0.0 {box_dim[0]}     xlo xhi\n")
            f.write(f"     -5.0 {box_dim[1]}     ylo yhi\n")
            f.write(f"     -5.0 {box_dim[2]}     zlo zhi\n\n")
            f.write("Masses\n\n")
            for i, t in enumerate(atom_types):
                f.write(f"     {i+1}   {self.type_atom_size(t)}\n")
            f.write('\n')
            f.write("Atoms\n\n")
            for a, atom in enumerate(atoms):
                f.write(f"     {a+1} 1 {atom[1]} {atom[0][0]:.4f} {atom[0][1]:.4f} {atom[0][2]:.4f}\n")
            f.write("\n")
            f.write("Bonds\n\n")
            [f.write(f"     {b+1} {bond[1]} {bond[0][0]} {bond[0][1]}\n")
                for b, bond in enumerate(bonds)]
            f.write("\n")
            f.write("Angles\n\n")
            [f.write(f"     {ang+1} {angle[1]} {angle[0][0]} {angle[0][1]} {angle[0][2]}\n")
                for ang, angle in enumerate(angles)]
            f.write("\n")
            f.write("Dihedrals\n\n")
            [f.write(f"     {d+1} {dih[1]} {dih[0][0]} {dih[0][1]} {dih[0][2]} {dih[0][3]}\n")
                for d, dih in enumerate(dihedrals)]
        print('[FINISHED WRITING]')


if __name__ == '__main__':
    tt = tools()
    tt.gen_atom_init()
