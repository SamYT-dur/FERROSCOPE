from pymatgen.symmetry import *
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.cif import CifWriter, CifParser
import re
import numpy as np
import os
from ccdc.io import EntryReader


from pyxtal import pyxtal
from gemmi import cif

def main(refcode):
    csd_reader = EntryReader('CSD') #find spacegroup
    entry = csd_reader.entry(refcode)
    spg = entry.crystal.spacegroup_symbol
    #bad_list = ['P31c', 'R3m', 'R3', 'R3c', 'R-3m'] #sometimes fail when reading into pymatgen: Some occupancies ([1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) sum to > 1! If they are within the occupancy_tolerance, they will be rescaled. The current occupancy_tolerance is set to: 1.0 Species occupancies sum to more than 1!
    try:
    #    if spg in bad_list: #fixes some of the issues by reading into pyxtal first. i think the problem is multiple atoms on one site
        pyx_struc = pyxtal()
        pyx_struc.from_seed('OutputCIFs/CSD/{}.cif'.format(refcode))
        low_struc = pyx_struc.to_pymatgen()
        #else: #this is how it should run without errors. All other space groups seem fine for now. not definitive
            #low_struc = Structure.from_file('OutputCIFs/CSD/{}.cif'.format(refcode), merge_tol=0.01) #can also compare to simplified cifs using 'MasterOutputs/Top/{}_top.cif' as string.
        high_struc = Structure.from_file('OutputCIFs/Output/{}_top_findsym.cif'.format(refcode))
        new_struc = Structure.from_sites(low_struc.sites)
    except Exception as error:
        with open('errors.txt', 'a+') as outfile:
            outfile.write(str(refcode) + ' : SUPERIMPOSE(cif parsing): ' + str(error) + '\n')
        return

    with open('OutputCIFs/Output/{}_top_findsym.cif'.format(refcode)) as infile:
        lines = infile.readlines()
        for i, line in enumerate(lines):
            if re.search('Origin', line):
                x_o_shift = float(line.split()[2])
                y_o_shift = float(line.split()[3])
                z_o_shift = float(line.split()[4])
            if re.search('Vectors a,b,c:', line):
                matrix_line = i

    trans_matrix = np.array([[float(lines[matrix_line+1].split()[0]), float(lines[matrix_line+1].split()[1]), float(lines[matrix_line+1].split()[2])],
                             [float(lines[matrix_line+2].split()[0]), float(lines[matrix_line+2].split()[1]), float(lines[matrix_line+2].split()[2])],
                             [float(lines[matrix_line+3].split()[0]), float(lines[matrix_line+3].split()[1]), float(lines[matrix_line+3].split()[2])]]
                            )
    origin_shift = np.array([x_o_shift, y_o_shift, z_o_shift])
    origin_shift = origin_shift.reshape(3, 1)
    trans_matrix2 = trans_matrix.transpose()

    high_species = []
    for species in high_struc.species:
        if species == 'C':
            high_species.append(species)
        else:
            high_species.append('O')
    high_coords2 = []

    for i in high_struc.frac_coords:
        high_coords2.append(i)
    high_coords3 = []

    for coords in high_coords2:
        coords2 = np.array([coords])
        coords2 = coords2.transpose()
        coords3 = np.matmul(trans_matrix2, coords2)
        coords4 = coords3 + origin_shift
        coords5 = [float(coords4[0]), float(coords4[1]), float(coords4[2])]
        high_coords3.append(coords5)

    for species, coords in zip(high_species, high_coords3):
        new_struc.append(species, coords, coords_are_cartesian=False)

    a = CifWriter(new_struc)
    a.write_file('OutputCIFs/Superimpose/{}_comp.cif'.format(refcode))


if __name__ == '__main__':
    with open('compcif_errors.txt', 'w+') as outfile:
        for file in os.listdir('OutputCIFs/Output'):
            filename = os.fsdecode(file)
            refcode = filename.strip('_top_findsym.cif')
           #try:
            main(refcode)
           #except:
            outfile.write(str(refcode) + '\n')
               #continue
