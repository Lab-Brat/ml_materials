import os
import numpy as np
from toolbox import tools
from pprint import pprint

class gen_chain():
    def __init__(self):
        # general settings
        self.t = tools()
        self.atom_num = 12
        self.init_atom = [5, 10, 10]
        self.box_dim = [50, 20, 20]

        # data dict. structure: 'object': [[obj, type], ...]
        self.polydat = {'atoms': [], 'bonds': [],
                        'angles': [], 'dihedrals': []}

        # define parameters for the algorithm
        self.elem_fun = {'C': 1}
        # self.elem_fun = {'C':0.5, 'N':0.5}
        self.elem_add = {'H': 1}
        self.elem_str = {'X': 1, 'Y': 2}

        # process parameters
        self.list_fun = [f for f in self.elem_fun.keys()]
        self.list_add = [a for a in self.elem_add.keys()]
        self.prob_fun = [self.elem_fun[i] for i in self.elem_fun]
        self.prob_add = [self.elem_add[i] for i in self.elem_add]
        self.types = {t: i+1 for i, t in enumerate(self.list_fun + self.list_add)}
        self.unit = self.elem_str['X']+self.elem_str['Y']
        self.blocks = self.atom_num // self.unit
        self.return_fun = lambda: np.random.choice(self.list_fun, p=self.prob_fun)
        self.return_add = lambda: np.random.choice(self.list_add, p=self.prob_add)

        # get split structure
        self.split = self.gen_str(self.atom_num)

    # polymer generating algorithm
    def gen_str(self, atom_num):
        structure = []
        for i in range(self.blocks):
            structure.append(self.return_fun())
            for a in range(self.elem_str['Y']):
                structure.append(self.return_add())
        split = [structure[i:self.unit+i] for i in range(0, atom_num, self.unit)]

        return split

    def mol_one(self):
        unit = self.split[0]
        bond1 = f'{unit[0]}{unit[1]}'
        bond2 = f'{unit[0]}{unit[2]}'
        self.atom = np.array(self.init_atom)
        H1, H2 = self.t.gen_atom_init([bond1, bond2], self.init_atom)
        at = [self.types[self.split[0][i]] for i in range(3)]

        self.polydat['atoms'].extend([[self.atom, at[0]],
                                      [H1, at[1]], [H2, at[2]]])
        self.polydat['bonds'].extend([[[1, 2], 2],
                                      [[1, 3], 2]])
        self.polydat['angles'].append([[2, 1, 3], 1])

    def _create_main_atom(self, mol):
        # C-C
        # print(i, mol)
        # bond_main = f'{mol[0]}{self.split[i][0]}'
        prev_main_ind = [i for i, e in enumerate(self.split) if e == mol][-1]
        bond_main = f'{mol[0]}{self.split[prev_main_ind][0]}'
        patom = self.polydat['atoms'][self.atid-3][0]
        atom = self.t.gen_atom_main(bond_main, patom)
        at = self.types[mol[0]]

        self.polydat['atoms'].extend([[atom, at]])
        self.polydat['bonds'].extend([ [[self.atid-2, self.atid+1], 1] ])

    def _create_side_atom(self, mol):
        # C-H
        last_c = self.polydat['atoms'][-1][0]
        bond1 = f'{mol[0]}{mol[1]}'
        bond2 = f'{mol[0]}{mol[2]}'
        H1, H2 = self.t.gen_atom_chain([bond1, bond2], last_c)
        at = [self.types[mol[i]] for i in range(1, 3)]

        self.polydat['atoms'].extend([[H1, at[0]],
                                      [H2, at[1]]])
        self.polydat['bonds'].extend([[[self.atid, self.atid+1], 2],
                                      [[self.atid, self.atid+2], 2]])
        self.polydat['angles'].append([[self.atid+1, self.atid, self.atid+2], 1])

    def mol_rest(self):
        for i, mol in enumerate(self.split[1:]):
            self.atid = len(self.polydat['atoms'])
            self.t.transpose()
            self._create_main_atom(mol)
            self.atid = len(self.polydat['atoms'])

            # add C-C-C angle and dihedrals
            if i > 0:
                self.polydat['angles'].append([[self.atid-6, self.atid-3, self.atid], 1])
            if i > 1:
                self.polydat['dihedrals'].append([[self.atid-9, self.atid-6,
                                                   self.atid-3, self.atid], 1])

            self._create_side_atom(mol)

    def build(self, write=False, copy=False):
        pprint(self.split)
        for i in range(1):
            self.mol_one()
            self.mol_rest()
            # self.init_atom = [10, 10, 10]

        if write is True:
            self.t.write_data(list(self.types.keys()),
                              self.polydat['atoms'],  self.polydat['bonds'],
                              self.polydat['angles'], self.polydat['dihedrals'],
                              self.box_dim)

        if copy is True:
            user = 'labbrat'
            home_path = f'/home/{user}'
            git_path = f'{home_path}/gitlab/polymer_gen'
            lmp_path = f'{home_path}/lammps/data/simple_minimization'
            os.system(f'cp {git_path}/carbon_chain.data \
                       {lmp_path}/polymerchain.data')
            print('[FINISHED COPYING]')

    def test(self):
        print(self.types[self.split[11][0]])


if __name__ == "__main__":
    gc = gen_chain()
    gc.build(write=True, copy=True)
    # gc.build(write=False, copy=False)
    # gc.test()
