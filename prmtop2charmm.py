import argparse
import re
import itertools
from math import degrees, pow
try:
    import networkx as nx
except:
    import sys
    sys.path.append('lib/')
    import networkx as nx

class molecules():
    title = None
    atom_list = []
    atoms = {}
    bonds = []
    angles = []
    torsions = []
    atom_type = {}
    params = {'bond': {'force': [], 'dist': []}, 'angle': {'force': [], 'phi': []}, 'torsion': {'force': [], 'period': [], 'phase': [], 'scee': [], 'scnb': []}, 'lj': {'a': [], 'b': []}}
    resnames = []
    residues = {}

class molecule(nx.Graph):
    _rtf = """
RESI %(resi)8s %(charge)9.4f
GROUP
%(atoms)s
%(bonds)s
"""
    _atom = "ATOM %(name)8s %(type)8s %(charge)9.4f"
    _bond = "BOND %s"

    def rtf(self):
        resi = self.graph['residue']
        atoms = []
        bonds = []
        charge = 0.0
        for atom, data in self.nodes(True):
            data['name'] = atom
            charge += data['charge']
            atoms.append(self._atom % data)

        edges = self.edges()
        for i in range(0, len(edges), 5):
            bonds.append(self._bond % " ".join(["%4s %4s" % edge for edge in edges[i:i+5]]))

        atoms = "\n".join(atoms)
        bonds = "\n".join(bonds)
        return self._rtf % locals()
        
class prmtop_parser():
    def parse(self, prmtopfile):
        fp = open(prmtopfile)
        flag = None
        self.mol = molecules()
        self.format = None
        for line in fp.readlines():
            line = line.rstrip()
            if line.startswith('%FLAG'):
                flag = line[6:]
                self.count = 0
                self.container = []
                continue
            if line.startswith('%FORMAT'):
                format = line[8:12]
                self.format = re.match("([0-9]+)([a-zA-Z])([0-9\.]+)", format).groups()
                continue

            if flag:
                f = getattr(self, 'parse_%s' % flag.lower())
                f(line)

    def parse_title(self, line):
        self.mol.title = line.strip()

    def parse_pointers(self, line):
        pass

    def parse_atom_name(self, line):
        for name in self._parse(line):
            self.mol.atom_list.append(name)
            self.mol.atoms[name] = {'id': self.count}
            self.count += 1

    def parse_charge(self, line):
        for charge in self._parse(line):
            name = self.mol.atom_list[self.count]
            self.mol.atoms[name]['charge'] = float(charge) / 18.2223
            self.count += 1

    def parse_atomic_number(self, line):
        pass

    def _parse(self, line):
        offset = int(self.format[2])
        n = int(self.format[0])
        for i in range(n):
            entry = line[i*offset:i*offset+offset].strip()
            if not entry: break
            yield entry

    def _parse_list(self, line, n):
        for entry in self._parse(line):
            self.container.append(entry)
            if len(self.container) == n:
                yield self.container
                self.container = []

    def parse_mass(self, line):
        for mass in self._parse(line):
            name = self.mol.atom_list[self.count]
            self.mol.atoms[name]['mass'] = float(mass)
            self.count += 1

    def parse_atom_type_index(self, line):
        for index in self._parse(line):
            name = self.mol.atom_list[self.count]
            self.mol.atoms[name]['type_index'] = int(index)
            self.count += 1

    def parse_excluded_atoms(self, line): pass
    def parse_number_excluded_atoms(self, line): pass
    def parse_nonbonded_parm_index(self, line): pass
    
    def parse_residue_label(self, line):
        for resn in self._parse(line):
            self.mol.residues[resn] = []
            self.mol.resnames.append(resn)

    def parse_residue_pointer(self, line):
        resid = [None, None]
        for index in self._parse(line):
            if resid[0]:
                resid[1] = int(index)
            else:
                resid[0] = int(index)
            if None not in resid:
                self._add_residue(resid)
                resid = [resid[1], None]

        self._add_residue((resid[0], len(self.mol.atom_list)+1))

    def _add_residue(self, resid):
        resname = self.mol.resnames[self.count]
        for i in range(resid[0]-1, resid[1]-1):
            self.mol.residues[resname].append(self.mol.atom_list[i])
            self.mol.atoms[self.mol.atom_list[i]]['resname'] = resname
            self.mol.atoms[self.mol.atom_list[i]]['resid'] = self.count + 1
        self.count += 1

    def parse_bond_force_constant(self, line):
        for k in self._parse(line):
            self.mol.params['bond']['force'].append(float(k))

    def parse_bond_equil_value(self, line):
        for dist in self._parse(line):
            self.mol.params['bond']['dist'].append(float(dist))

    def parse_angle_force_constant(self, line):
        for k in self._parse(line):
            self.mol.params['angle']['force'].append(float(k))

    def parse_angle_equil_value(self, line):
        for phi in self._parse(line):
            self.mol.params['angle']['phi'].append(degrees(float(phi)))

    def parse_dihedral_force_constant(self, line):
        for k in self._parse(line):
            self.mol.params['torsion']['force'].append(float(k))

    def parse_dihedral_periodicity(self, line):
        for period in self._parse(line):
            self.mol.params['torsion']['period'].append(float(period))

    def parse_dihedral_phase(self, line):
        for phi in self._parse(line):
            self.mol.params['torsion']['phase'].append(degrees(float(phi)))

    def parse_scee_scale_factor(self, line):
        for scale in self._parse(line):
            self.mol.params['torsion']['scee'].append(float(scale))

    def parse_scnb_scale_factor(self, line):
        for scale in self._parse(line):
            self.mol.params['torsion']['scnb'].append(float(scale))

    def parse_lennard_jones_acoef(self, line):
        for sigma in self._parse(line):
            self.mol.params['lj']['a'].append(float(sigma))

    def parse_lennard_jones_bcoef(self, line):
        for sigma in self._parse(line):
            self.mol.params['lj']['b'].append(float(sigma))

    def _convert_atom_list(self, l):
        return [int(x)/3 if i<(len(l)-1) else int(x) for i,x in enumerate(l)]

    def parse_bonds_inc_hydrogen(self, line):
        for bond in self._parse_list(line, 3):
            self.mol.bonds.append(self._convert_atom_list(bond))

    def parse_bonds_without_hydrogen(self, line):
        for bond in self._parse_list(line, 3):
            self.mol.bonds.append(self._convert_atom_list(bond))

    def parse_angles_inc_hydrogen(self, line):
        for angl in self._parse_list(line, 4):
            self.mol.angles.append(self._convert_atom_list(angl))

    def parse_angles_without_hydrogen(self, line):
        for angl in self._parse_list(line, 4):
            self.mol.angles.append(self._convert_atom_list(angl))

    def parse_dihedrals_inc_hydrogen(self, line):
        for tors in self._parse_list(line, 5):
            self.mol.torsions.append(self._convert_atom_list(tors))

    def parse_dihedrals_without_hydrogen(self, line):
        for tors in self._parse_list(line, 5):
            self.mol.torsions.append(self._convert_atom_list(tors))

    def parse_amber_atom_type(self, line):
        for atomtype in self._parse(line):
            name = self.mol.atom_list[self.count]
            self.count += 1
            index = self.mol.atoms[name]['type_index']
            self.mol.atoms[name]['type'] = atomtype
            if self.mol.atom_type.has_key(index): continue
            self.mol.atom_type[index] = {'id': index, 'type': atomtype, 'mass': self.mol.atoms[name]['mass']}

    def parse_solty(self, line): pass
    def parse_hbond_acoef(self, line): pass
    def parse_hbond_bcoef(self, line): pass
    def parse_hbcut(self, line): pass
    def parse_excluded_atoms_list(self, line): pass
    def parse_tree_chain_classification(self, line): pass
    def parse_join_array(self, line): pass
    def parse_irotat(self, line): pass
    def parse_ipol(self, line): pass
    def parse_screen(self, line): pass
    def parse_radii(self, line): pass
    def parse_radius_set(self, line): pass

def build_molecules(prmtop):
    molecules = []
    prms = {'bond': {}, 'angle': {}, 'torsion': {}, 'lj': {}}
    for i, resname in enumerate(set(prmtop.resnames)):
        resi = molecule(residue=resname)
        for atom in prmtop.residues[resname]:
            data = prmtop.atoms[atom]
            resi.add_node(atom, data)

        for bond in prmtop.bonds:
            atom1 = prmtop.atom_list[bond[0]]
            atom2 = prmtop.atom_list[bond[1]]
            if resi.has_node(atom1) and resi.has_node(atom2):
                resi.add_edge(atom1, atom2)
        molecules.append(resi)

    for bond in prmtop.bonds:
        atom1 = prmtop.atom_list[bond[0]]
        atom2 = prmtop.atom_list[bond[1]]

        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        pair = tuple(sorted((type1, type2)))
        if prms['bond'].has_key(pair): continue
        index = bond[2] - 1
        prms['bond'][pair] = \
            {'force': prmtop.params['bond']['force'][index],
             'dist': prmtop.params['bond']['dist'][index]}

    for angl in prmtop.angles:
        atom1 = prmtop.atom_list[angl[0]]
        atom2 = prmtop.atom_list[angl[1]]
        atom3 = prmtop.atom_list[angl[2]]

        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        type3 = prmtop.atoms[atom3]['type']
        pair = tuple(sorted((type1, type2, type3)))
        if prms['bond'].has_key(pair): continue
        index = angl[3] - 1
        prms['angle'][pair] = \
            {'force': prmtop.params['angle']['force'][index],
             'phi': prmtop.params['angle']['phi'][index]}

    for tors in prmtop.torsions:
        atom1 = prmtop.atom_list[tors[0]]
        atom2 = prmtop.atom_list[tors[1]]
        atom3 = prmtop.atom_list[tors[2]]
        atom4 = prmtop.atom_list[tors[3]]

        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        type3 = prmtop.atoms[atom3]['type']
        type4 = prmtop.atoms[atom4]['type']
        pair = tuple(sorted((type1, type2, type3, type4)))
        if prms['torsion'].has_key(pair): continue
        index = tors[4] - 1
        prms['torsion'][pair] = \
            {'force': prmtop.params['torsion']['force'][index],
             'phase': prmtop.params['torsion']['phase'][index],
             'period': prmtop.params['torsion']['period'][index],
             'scee': prmtop.params['torsion']['scee'][index],
             'scnb': prmtop.params['torsion']['scnb'][index]}

    for i, pair in enumerate(itertools.combinations(sorted(prmtop.atom_type.keys()), 2)):
        prms['lj'][pair] = {'a': prmtop.params['lj']['a'][i], 'b': prmtop.params['lj']['b'][i]}

    return molecules, prms

def build_rtf(molecules, prmtop, prefix):
    _rtf = """* %(title)s
*
36 1

%(mass)s

%(residues)s
"""
    _mass = "MASS %(type_index)8d %(type)4s %(mass)9.5f"
    atoms = []

    residues = ""
    for mol in molecules:
        atoms.extend([_mass % data for atom,data in mol.nodes(True)])
        residues += mol.rtf()

    mass = "\n".join(sorted(set(atoms)))
    title = prmtop.title

    fp = open('%s.rtf' % prefix, 'w')
    fp.write(_rtf % locals())

def build_prm(prmtop, prefix):
    _prm = """* %(title)s
*

ATOM
%(mass)s

BONDS
%(bonds)s

ANGLES
%(angles)s

DIHEDRALS
%(torsions)s

NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 

NBFIX
%(nonbonds)s
"""
    _mass = "MASS %(type_index)8d %(type)4s %(mass)9.5f"
    _bond = "%(atom1)4s %(atom2)4s %(force)8.4f %(dist)8.4f"
    _angl = "%(atom1)4s %(atom2)4s %(atom3)4s %(force)8.4f %(phi)8.4f"
    _dihe = "%(atom1)4s %(atom2)4s %(atom3)4s %(atom4)4s %(force)8.4f %(period)4d %(phase)8.4f"
    _nonb = "%(type1)4s %(type2)4s %(eps)16.8f %(rmin)16.8f"

    atoms = [_mass % {'type_index': k, 'type': v['type'], 'mass': v['mass']} for k,v in prmtop.atom_type.items()]

    bonds = {}
    for prm in prmtop.bonds:
        if bonds.has_key(prm[2]): continue
        atom1 = prmtop.atom_list[prm[0]]
        atom2 = prmtop.atom_list[prm[1]]
        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        bonds[prm[2]] = {'atom1': type1, 'atom2': type2, 'force': prmtop.params['bond']['force'][prm[2]-1], 'dist': prmtop.params['bond']['dist'][prm[2]-1]}

    angles = {}
    for prm in prmtop.angles:
        if angles.has_key(prm[3]): continue
        atom1 = prmtop.atom_list[prm[0]]
        atom2 = prmtop.atom_list[prm[1]]
        atom3 = prmtop.atom_list[prm[2]]
        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        type3 = prmtop.atoms[atom3]['type']
        angles[prm[3]] = {'atom1': type1, 'atom2': type2, 'atom3': type3, 'force': prmtop.params['angle']['force'][prm[3]-1], 'phi': prmtop.params['angle']['phi'][prm[3]-1]}

    torsions = {}
    for prm in prmtop.torsions:
        #if torsions.has_key(prm[4]): continue
        atom1 = prmtop.atom_list[prm[0]]
        atom2 = prmtop.atom_list[prm[1]]
        atom3 = prmtop.atom_list[prm[2]]
        atom4 = prmtop.atom_list[prm[3]]
        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        type3 = prmtop.atoms[atom3]['type']
        type4 = prmtop.atoms[atom4]['type']
        torsions[prm[4]] = {'atom1': type1, 'atom2': type2, 'atom3': type3, 'atom4': type4, 'force': prmtop.params['torsion']['force'][prm[4]-1], 'phase': prmtop.params['torsion']['phase'][prm[4]-1], 'period': prmtop.params['torsion']['period'][prm[4]-1]}

    vdws = []
    for k, (i, j) in enumerate(itertools.combinations(prmtop.atom_type, 2)):
        atom1 = prmtop.atom_list[i]
        atom2 = prmtop.atom_list[j]
        type1 = prmtop.atoms[atom1]['type']
        type2 = prmtop.atoms[atom2]['type']
        a = prmtop.params['lj']['a'][k]
        b = prmtop.params['lj']['b'][k]
        sigma = pow(a/b, 1/6.0)
        eps = -b**2 / 4.0 / a
        rmin = pow(2, 1/6.0)*sigma
        vdws.append(_nonb % locals())

    mass = "\n".join(sorted(set(atoms)))
    bonds = "\n".join([_bond % bonds[k] for k in sorted(bonds.keys())])
    angles = "\n".join([_angl % angles[k] for k in sorted(angles.keys())])
    torsions = "\n".join([_dihe % torsions[k] for k in sorted(torsions.keys())])
    nonbonds = "\n".join(vdws)
    title = prmtop.title

    fp = open('%s.prm' % prefix, 'w')
    fp.write(_prm % locals())

def main():
    parser = argparse.ArgumentParser(description='translate prmtop file to psf, rtf, and prm')    
    parser.add_argument('prmtop', help='prmtop filename')
    parser.add_argument('prefix', help='prefix name for output')
    args = parser.parse_args()

    parser = prmtop_parser()
    parser.parse(args.prmtop)
    molecules, prms = build_molecules(parser.mol)

    build_rtf(molecules, parser.mol, args.prefix)
    build_prm(parser.mol, args.prefix)

if __name__ == '__main__':
    main()
