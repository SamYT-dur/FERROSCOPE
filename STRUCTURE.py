# script for topological labelling of molecules to find equivalent atoms

import sys
import math
import os
import glob
import stopit

from ccdc.io import CrystalReader
from ccdc.io import CrystalWriter
from ccdc.io import EntryReader
from ccdc.molecule import Bond
import TOOLS
import ccdc.io, ccdc.search, ccdc.entry

# File set-up for running standalone. Otherwise, these are set by MASTER.py
# region
replace_type_by_c = ['N', 'O']  # atoms that get changed to C because they'll be similar in size and could disorder
excluded_atoms = ['H', 'D']  # atoms to delete before symmetry analysis etc
ABn_list = ['C_Cl_3', 'C_F_3', 'C_Br_3', 'Cl_O_4', 'Br_O_4', 'I_O_4', 'S_O_4', 'Cl_C_4', 'Br_C_4', 'I_C_4', 'B_F_4', 'B_H_4', 'Fe_Cl_4', 'Fe_Cl_5', 'Fe_Br_4', 'Re_O_4', 'Re_C_4', 'Ga_Cl_4', 'Mn_Br_4', 'Mn_Cl_4', 'Zn_Cl_4', 'Cu_Cl_4', 'Co_Cl_4', 'Al_Cl_4', 'Tl_Cl_4', 'Be_Cl_4', 'Cd_Cl_4', 'Zn_Br_4', 'Co_Br_4', 'Al_Br_4', 'Hg_Br_4', 'Ga_Br_4', 'Zn_I_4', 'Hg_I_4', 'In_I_4', 'Al_I_4']  # groups that likely to be disordered so B atoms removed - Sam -expanded by me, there could be more
notes_file = 'Outputs/process_notes.txt'
topology_relabel = True
# endregion


def bond_type_as_index(type):
    if type == Bond.BondType(1):  # single bond
        return 1
    if type == Bond.BondType(2):  # double
        return 2
    if type == Bond.BondType(3):  # triple
        return 3
    if type == Bond.BondType(4):  # quadruple
        return 4
    if type == Bond.BondType(5):  # aromatic
        return 5
    if type == Bond.BondType(7):  # delocalised
        return 7
    return 9  # pi bond


def calculate_next_level(molecule,current_level, strrep):
    """This fn calculates the equivalent of A(i, j) for each atom"""
    level = [ 0 for atom in molecule.atoms ]
    for atom in molecule.atoms:
        sum = 0
        for bond in atom.bonds:
            for batom in bond.atoms:
                if batom.index != atom.index:  # checks that the atom is not itself in the bond
                     #print(current_level)
                     #print(batom.index)
                     sum = sum + current_level[batom.index]


        strrep[atom.index] = str(strrep[atom.index]) + "_" + str(sum)
        level[atom.index] = sum
    return level,strrep


def calculate_initial_indexes(molecule):
    """ I think this fn just provides a way of separating out different atoms/bonding
so that they do not get assigned to the same groups during the partitions"""
    level = [ 0 for atom in molecule.atoms ]
    strrep = [ "" for atom in molecule.atoms ]
    for atom in molecule.atoms:
        sum = atom.atomic_number
        for bond in atom.bonds:
            sum = sum
                  # + math.pow(10,bond_type_as_index(bond.bond_type) + 2)  # comment this out so we ignore bonding in deciding equivalence
        strrep[atom.index] = str(sum)
        level[atom.index] = sum
    return level,strrep


def calculate_indexes(molecule):

    group_by_index = {}

    level,strrep = calculate_initial_indexes(molecule)

    n_partitions = len(set(strrep))  # finds the number of unique elements in strrep, Ng in the paper

    next_level,strrep = calculate_next_level(molecule,level,strrep)

    while len(set(strrep)) > n_partitions:  # iterates until the number of groups stops increasing
        n_partitions = len(set(strrep))
        next_level,strrep = calculate_next_level(molecule,next_level,strrep)
    return next_level,strrep


def output_crystal(identifier):
    csd_reader = EntryReader('csd')
    entry_id = identifier
    entry = csd_reader.entry(entry_id)
    TOOLS.entry_id = entry_id
    crystal = entry.crystal

    disordered = crystal.disordered_molecule  # first we take the full molecule, with disorder

    suppressed = [atom for atom in disordered.atoms if atom.label.endswith('?')]  # create a list of suppressed atoms

    new_mol = crystal.asymmetric_unit_molecule

    for atom in suppressed:  # here we remove any suppressed atoms from the asymmetric unit
        label = atom.label
        new_mol.remove_atom(new_mol.atom(label))

    molecule = new_mol  # rename new_mol to fit with the rest of the code
    note = ''
    num_simplified = 0

    # Options decided in MASTER.py
    if excluded_check:
        molecule = TOOLS.exclude(molecule,excluded,replace_type_by_c)  # exclude and replace atoms

    if solvents:
        molecule, note, num_simplified = TOOLS.solvents(molecule, note, num_simplified)  # Sam's second attempt at simplification of small organics using moment of inertia

    if tds:
        molecule, note = TOOLS.sam_Abn_A(molecule, crystal, ABn_list, note)  # sam's replacement for replace_ABn_A #operating on the asym unit does miss a couple entries but not worrying about that for now

    if len(note) > 0:
        print('**Output note for Results2.csv for entry',entry_id,'is:',note)
        with open(notes_file, 'a') as outfile:  # write the note to file
            outfile.write('{0}:{1}\n'.format(entry_id,note))
    if len(note) == 0:
        note='None'

# this bit here does the topology atom relabelling

    if topology_relabel:
        indexes, strrep = calculate_indexes(molecule)
        lookup = []
        for atom in molecule.atoms:
            lookup.append([strrep[atom.index],atom.index])
        lookup.sort()  # creates a 2D list of strreps and the corresponding atom number, then sorts them highest first
        group_counter = 0
        last_id = ""
        new_atoms = [None for atom in molecule.atoms]
        for ordered in lookup:
            atom = molecule.atoms[ordered[1]]  # picks out the first atom's index (ordered is a tuple)
            if ordered[0] != last_id:  # checks it is not the same as the last atom checked
                last_id = ordered[0]
                group_counter = group_counter + 1  # as it is different, we have a new unique group
            atom.label = atom.label + "_" + str(group_counter)  # relabels the atom with its group number, so atoms in the same group are related.
            new_atoms[ordered[1]] = atom

    crystal.molecule = molecule
    return crystal

if __name__ == '__main__':
    #id = input('Enter an entry ID from the CSD:')
#    id = 'FIGNER'
    open(notes_file, 'w').close() #empty the process notes file only needed when running this
    #print('****just deleted notes_file in STRUCTURE.py')
    idlist=['AACANI10','UBOTIQ01']
    csd_reader = ccdc.io.EntryReader('CSD')
    for entry in csd_reader:
        crystal = entry.crystal
        id = crystal.identifier
    #for id in idlist:
        #output_crystal(identifier=id)

        altered = output_crystal(identifier=id)
        top_an_cif = altered.to_string(format='cif')
        #print(top_an_cif)
        altered_out = id+'_top.cif'
    #    altered_out = os.path.join(top_an_out_dir, altered_file_name)
        with open(altered_out, 'w+') as outfile:
            outfile.write(top_an_cif)


        top_an_cif = altered.to_string(format='cif')
        #print('done identifier',id)
