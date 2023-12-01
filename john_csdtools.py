#
# some csd routines
# these can be called from e.g. john_top.py or EDG_topological_analysis.py
# started 15/5/2020
# 23/5/2020 added exclude CF3 type group via a C_F_3 type format in a list
# 27/1/2022 added new function simplify_small_organics
# 19/5/2023 sam testing replace_ABn_A. replaced this with sam_ABn_A

import math
import os
from iodata import load_one, dump_one
import ccdc.molecule
import ccdc.search
import ccdc.io
import ccdc.conformer
from ccdc.descriptors import MolecularDescriptors as MD
from ccdc.descriptors import GeometricDescriptors as GD

import traceback
import sys
import numpy as np
import sqlite3


#Imports
#------------------------------
#External Libraries
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.analyzer import generate_full_symmops

#PyXtal imports
from pyxtal.operations import *
from pyxtal.molecule import *
from pyxtal.viz import *
from pyxtal.database.collection import Collection

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D
from rdkit import DataStructs


def addtwo(waters):
    waters +=2
    return waters

def exclude(molecule,excluded,replace_type_by_c):
    #excluded = ['H','Ni','D','Cl','Br','O'] #normally set by EDG_Master; included here so EDG_Basic and this script will run standalone
    #replace_type_by_c = ['N','O','U'] #normally set by EDG_Master; included here so EDG_BAsic and this script will run standalone
    for atom in molecule.atoms:
        if atom.atomic_symbol in replace_type_by_c:
            atom.atomic_symbol = 'C'
        if atom.atomic_symbol in excluded:
            molecule.remove_atom(atom)
    return molecule


def remove_water(molecule,waters,note):
    keep = []
    waters = 0
    #print('entered remove_water function')
    #print('number of components is',len(molecule.components))
    for s in molecule.components:
        ats = [at.atomic_symbol for at in s.atoms]
        if len(ats) == 3:
            ats.sort()
            if ats[0] == 'H' and ats[1] == 'H' and ats[2] == 'O':
                waters += 1
                #print('found a water',ats[0],ats[1],ats[2])
            else:
                keep.append(s)
        elif len(ats) == 1: #sam added this for waters with no hydrogens
            if ats[0] == 'O':
                waters += 1
            else:
                keep.append(s)
        else:
            keep.append(s)
    new = ccdc.molecule.Molecule(molecule.identifier)
    for k in keep:
        new.add_molecule(k)
    #print('waters found',waters)
    molecule=new
    #print('number of components after waters removed',len(molecule.components))
    if waters > 0:
        note=note+'Excluded '+str(waters)+' waters; '
    return molecule,waters,note

#looks for AB4 or A4B species
#ought to put some check to see if tetrahedral
#this only works if the group in question has been correctly identified as a separate component in the database
def remove_xy4(molecule,xy4_exclude,tet_anion,note):
    keep = []
    tet_anion = 0
    #print('entered remove_tet_anion function')
    #print('number of components is',len(molecule.components))
    # - = ['ClO4','BF4','IO4']
    #note=''
    for s in molecule.components:
        ats = [at.atomic_symbol for at in s.atoms]
        if len(ats) == 5:
            ats.sort()
            if (ats[0] ==  ats[1] ==  ats[2] == ats[3]) or (ats[1] ==  ats[2] ==  ats[3] == ats[4]):
                #print('found a ab4 or a4b',ats[0],ats[1],ats[2],ats[3],ats[4])
                if ats[0] == ats[1]:
                    group_name=ats[1]+ats[0]+'4'
                else:
                    group_name=ats[0]+ats[1]+'4'
                #print('AB4 or A4B group found:',group_name)
                if group_name in xy4_exclude:
                    #print('I should exclude',group_name)
                    #print(note)
                    tet_anion += 1
                    #could put an if group_name in group_name_list to remove ClO4 IO4 BH4 etc and write to "note" for reporting
                    for atom in s.atoms:
                        if atom.atomic_symbol == ats[1]:
                            s.remove_atom(atom)
                            #print('removing',atom)
                keep.append(s)
        else:
            keep.append(s)
    new = ccdc.molecule.Molecule(molecule.identifier)
    if tet_anion >0:
        note=note+'Excluded '+str(tet_anion)+' AB4 probably all '+group_name+';'
    for k in keep:
        new.add_molecule(k)
    #print('ab4 or a4b found',tet_anion)
    molecule=new
    #print('number of components after ab4 removed',len(molecule.components))
    return molecule,tet_anion,note

#looks for AB4 or A4B species
#version 2 goes to each atom and looks through neighbours as some ClO4 groups not in their own component in database
#this is hugely inefficient as you end up looking at lots of irrelevant atoms
#10/6/2020 should be deleted from toolkit
def remove_xy4_v2(molecule,xy4_exclude,tet_anion,note):
    keep = []
    tet_anion = 0
    group_name = ''
    for s in molecule.components:
        for at in s.atoms:
            if len(at.neighbours) == 4:
                ats = [at.atomic_symbol for at in at.neighbours]
                if ats[0] == ats[1] == ats[2] == ats[3]:
                    group_name = at.atomic_symbol+ats[0]+'4'
                    if group_name in xy4_exclude:
                        tet_anion += 1
                        for at in at.neighbours:
                            s.remove_atom(at)
                            #print('deleting:',at)
            keep.append(s)
        else:
            keep.append(s)
    new = ccdc.molecule.Molecule(molecule.identifier)
    if tet_anion > 0:
        note=note+'Excluded '+str(tet_anion)+' AB4 probably all '+group_name+';'
    for k in keep:
        new.add_molecule(k)
    molecule=new
    return molecule,tet_anion,note



#looks for species A_B_n like C_Cl_3 in a molecule and deletes the B atoms
#ABn_list=['C_Cl_3','N_H_2']
def delete_ABn(molecule, ABn_list, note):
    for group in ABn_list:
        ABn=0
        keep=[]
        centre,bonded,cnum=group.split('_')
        #print(centre,bonded,cnum)
        for s in molecule.components:
            for at in s.atoms:
                neighlist = [at.atomic_symbol for at in at.neighbours]
                if at.atomic_symbol == centre and neighlist.count(bonded) == int(cnum):
                    ABn +=1
                    for neigh in at.neighbours:
                        if neigh.atomic_symbol == bonded:
                            s.remove_atom(neigh)
            keep.append(s)
        new = ccdc.molecule.Molecule(molecule.identifier)
        for k in keep:
            new.add_molecule(k)
        molecule=new
        if ABn > 0:
            note=note+'Removed '+str(ABn)+'*'+group.replace('_','')+'; '
    return molecule, note

#
# Use CSD search tools to find ABn species which should be more robust
# Will search for ABn bonded group (any bonds) and remove all the B
#
def replace_ABn_A(molecule, ABn_list, note):
    match_labels = []
    keep=[]
    for group in ABn_list:
        #print('working with group',group)
        centre,bonded,cnum=group.split('_')
        #create a substructure from Cl_O_4 name
        ABn_substructure = ccdc.search.QuerySubstructure()
        exec('a1 = ABn_substructure.add_atom(centre)')

        for x in range(1,int(cnum)+1):
            exec('b'+str(x)+' = ABn_substructure.add_atom("'+bonded+'")')
            exec('bo'+str(x)+' = ABn_substructure.add_bond("Any", a1, b'+str(x)+')')
        ABn_search = ccdc.search.SubstructureSearch()
        #ABn_search.add_atom_property_constraint('ATOM1', (0,0), int(cnum), which='neighbours')
        sub_id= ABn_search.add_substructure(ABn_substructure)

        hits = ABn_search.search([molecule.identifier])
        if len(hits) > 0:
            note=note+'Removed '+str(len(hits))+'*'+group.replace('_','')+'; '

        #get all search matches into a single list of atoms to remove which dont contain central atom
        for hit in hits:
            matches=hit.match_atoms()
            print('atom matches:',matches)
            match_labels = match_labels+[match.label for match in matches]
        #this gets all neighbours that contain the bonded atom.  Danger that O would remove O1 and Os1 so use a re.sub
        match_labels[:] = [x for x in match_labels if centre not in x] #This is broken. CCl3 doesnt work because C is in Cl

    #print('full list of atomic labels to exclude from all groups is',match_labels)

    for s in molecule.components:
        for at in s.atoms:
            if at.label in match_labels:
                s.remove_atom(at)
        keep.append(s)
    new = ccdc.molecule.Molecule(molecule.identifier)
    for k in keep:
        new.add_molecule(k)
    molecule = new
    return molecule, note


#Sam's attempt at replacing replace_ABn_a. fixes a lot of issues with the previous iteration but has problems of it's own. Operating on the asymmetric unit means that some tetrahedrons will look like A-B2 rather than A-B4 for example and won't be simplified.
def sam_Abn_A(molecule, crystal, ABn_list, note):
    to_remove = []
    crystmol = crystal.molecule
    for group in ABn_list:
        centre, bonded, cnum = group.split('_')
        for comp in crystmol.components:
            for atom in comp.atoms:
                if atom.atomic_symbol == centre and atom.atomic_symbol != 'C':
                    bondeds = []
                    for a in atom.neighbours:
                        if a.atomic_symbol == bonded and len(atom.bonds) == int(cnum):
                            bondeds.append(a.label)
                    if len(bondeds) == int(cnum):
                        for l in bondeds:
                            to_remove.append(l)
                if atom.atomic_symbol == centre and atom.atomic_symbol == 'C':
                    bondeds = []
                    for a in atom.neighbours:
                        if a.atomic_symbol == bonded and len(atom.bonds) == int(cnum) + 1:
                            bondeds.append(a.label)
                    if len(bondeds) == int(cnum):
                        for l in bondeds:
                            to_remove.append(l)
    for atom in molecule.atoms:
        if atom.label in to_remove:
            molecule.remove_atom(atom)

    return molecule, note

#Sam - straightforward removal of any molecule containing a DABCO-like substructure.
def simplify_dabco(molecule, note, num_simplified):

    dabco = Chem.MolFromSmiles('C1CC2CCC1CC2')
    for component in molecule.components:
        if len(component.atoms) < 15:
            smiles = component.smiles
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                overlap = mol.GetSubstructMatch(dabco)
                if overlap != ():
                    centroid_coords = MD.atom_centroid(*tuple(a for a in component.atoms))
                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]), label='He100{}'.format(str(num_simplified)))
                    new_component.add_atom(atom)
                    for comp2 in molecule.components:
                        if comp2 != component:
                            new_mol.add_molecule(comp2)
                    new_mol.add_molecule(new_component)
                    molecule = new_mol
                    num_simplified += 1

            except:
                pass
    return molecule, note, num_simplified

#Sam - this function is now made redundant by inertial below
def simplify_small_organics(molecule, note, num_simplified):
    if molecule.is_organometallic:
        for component in molecule.components:
            if len(component.atoms) < 14:
                carbons = 0
                other_atoms = 0
                for atom in component.atoms:
                    atom_remove = ['F', 'Cl', 'Br', 'I']
                    if atom.atomic_symbol == 'C':
                        carbons += 1
                    if atom.atomic_symbol in atom_remove:
                        other_atoms += 1
                if carbons + other_atoms == len(component.atoms) and other_atoms != len(component.atoms):  #added to prevent single chloride ions being turned into helium
                    centroid_coords = str(MD.atom_centroid(*tuple(a for a in component.atoms)))
                    x_place = centroid_coords.find('x')
                    centroid_x = centroid_coords[x_place+2:x_place+9]
                    centroid_x = centroid_x.replace(',','')
                    y_place = centroid_coords.find('y')
                    centroid_y = centroid_coords[y_place+2:y_place+9]
                    centroid_y = centroid_y.replace(',','')
                    z_place = centroid_coords.find('z')
                    centroid_z = centroid_coords[z_place+2:z_place+9]
                    centroid_z = centroid_z.replace(',','')
                    centroid_z = centroid_z.replace(')', '')
                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                    new_mol = ccdc.molecule.Molecule(identifier='new_mol')
                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(float(centroid_x), float(centroid_y), float(centroid_z)), label='He100')
                    new_component.add_atom(atom)
                    for comp in molecule.components:
                        if comp != component:
                            new_mol.add_molecule(comp)
                    new_mol.add_molecule(new_component)
                    molecule = new_mol
    return (molecule, note, num_simplified)


def simplify_small_organics_inertial(molecule, note, num_simplified): #Sam's V2
        for component in molecule.components:
                #print(component.formula)
                #print('xoxoxoxoxoxoxoxox' + component.identifier)
                if len(component.atoms) < 14: #should be 14 for normal runs - Sam
                    carbons = 0
                    other_atoms = 0
                    for atom in component.atoms:
                        atom_remove = ['F', 'Cl', 'Br', 'I']
                        if atom.atomic_symbol == 'C':
                            carbons += 1
                        if atom.atomic_symbol in atom_remove:
                            other_atoms += 1
                    if carbons + other_atoms == len(component.atoms) and other_atoms != len(component.atoms):
                        small_organics_reader = ccdc.io.CrystalReader('CSD')  #the next lines pick out the pure molcule from the database to avoid asymmetric unit only being part of a molecule. Also undoes any simplification to enable more accurate inertia values.
                        small_organics_molecule = small_organics_reader.molecule('{}'.format(molecule.identifier))
                        for comp in small_organics_molecule.components:
                            if comp.identifier == component.identifier:
                                try:  # incase adding hydrogens fails e.g. unknown bond type
                                    comp.add_hydrogens(mode='missing', add_sites=True)
                                    #molecule_minimiser = ccdc.conformer.MoleculeMinimiser() ##currently not minimising so inertia values are non-minimised to accurately represent conformer in crystal
                                    #comp2 = molecule_minimiser.minimise(comp)
                                    #comp = comp2
                                except RuntimeError:
                                    print('failed to add hydrogens')
                                    for atom in comp.atoms: # prevents 'atom has no coordinates error'
                                        if atom.coordinates == None:
                                            comp.remove_atom(atom)
                                    pass
                                mol_file = comp.to_string(format='mol2')  #csd doesn't output .xyz files so output and write out mol2 file
                                if os.path.isfile('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier))) == False:
                                    with open('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)), 'w+') as datafile:
                                        datafile.write(mol_file)
                                else:

                                    with open('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)), 'w+') as datafile:
                                        datafile.write(mol_file)
                                try:
                                    iodata_mol = load_one('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)))  #these two lines us iodata module to convert mol2 to xyz and write out
                                    dump_one(iodata_mol, 'small_mols/xyz/{}{}.xyz'.format(small_organics_molecule.identifier, str(comp.identifier)))
                                    pymatgen_mol = mol_from_file('small_mols/xyz/{}{}.xyz'.format(small_organics_molecule.identifier, str(comp.identifier)))
                                    reo_mol = reoriented_molecule(pymatgen_mol)[0]
                                    rotation_matrix = reoriented_molecule(pymatgen_mol)[1].rotation_matrix
                                    reo_mol.to(filename='small_mols/reo/{}.json'.format(small_organics_molecule.identifier))
                                    x_inertia = get_moment_of_inertia(reo_mol, [1, 0, 0])
                                    y_inertia = get_moment_of_inertia(reo_mol, [0, 1, 0])
                                    z_inertia = get_moment_of_inertia(reo_mol, [0, 0, 1])
                                    spherocity_index = (3*x_inertia)/(x_inertia+y_inertia+z_inertia) #with weights
                                except Exception as e:
                                    with open('pymat_errs.txt', 'a+') as out:
                                        out.write(str(molecule.identifier) + ': pymatgen error --> {}\n'.format(e))
                                    #print('pymatgen error')
                                    break
                                try:
                                    pyx_mol = pyxtal_molecule(reo_mol)
                                    box = pyx_mol.get_box(padding=0)
                                    box_length = box.length
                                except ValueError as err: #converting to pyx mol sometimes fails for ammonium
                                    #print(traceback.format_exc())
                                    atoms = []
                                    for atom in comp.atoms:
                                        atoms.append(atom.atomic_symbol)
                                        if atoms == ['N', 'H', 'H', 'H', 'H']:
                                            new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                            for comp2 in molecule.components:
                                                if comp2 != component:
                                                    new_mol.add_molecule(comp2)
                                            molecule = new_mol
                                            note = note + 'Removed {}; '.format(comp.formula)
                                            break
                                    with open('pymat_errs.txt', 'a+') as out:
                                        out.write(str(molecule.identifier) + ': pyxtal error\n')
                                    #print('pxtal error')
                                    break

                                #with open('spherocity indecies.txt', 'a+') as out:
                                #    out.write(molecule.identifier + ':  ' + str(spherocity_index) + '\n')

                                if x_inertia < 3 and y_inertia < 3 and z_inertia < 5: # removing water and ammonium
                                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                    for comp2 in molecule.components:
                                        if comp2 != component:
                                            new_mol.add_molecule(comp2)
                                    molecule = new_mol
                                    note = note + 'Removed {}; '.format(comp.formula)
                                    break

                                if 0.94 <= spherocity_index <= 1.06:   #simplify spherical molecules by replacing with a helium at the centroid
                                    try:
                                        centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                        new_component = ccdc.molecule.Molecule(identifier='centroid')
                                        new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                        atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]), label='He100{}'.format(str(num_simplified)))
                                        new_component.add_atom(atom)
                                        for comp2 in molecule.components:
                                            if comp2 != component:
                                                new_mol.add_molecule(comp2)
                                        new_mol.add_molecule(new_component)
                                        note = note + 'Replaced spherical molecule {} with point; '.format(comp.formula)
                                        molecule = new_mol
                                        num_simplified += 1
                                        break
                                    except:
                                        pass

                                #print(component.formula)
                                if (x_inertia <= 45 or y_inertia <= 45) and molecule.is_organometallic == True:
                                    centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]),label='He100{}'.format(str(num_simplified)))
                                    new_component.add_atom(atom)
                                    for comp2 in molecule.components:
                                        if comp2 != component:
                                            new_mol.add_molecule(comp2)
                                    new_mol.add_molecule(new_component)
                                    molecule = new_mol
                                    num_simplified += 1
                                    break
                                if (x_inertia <= 25 or y_inertia <= 25) and molecule.is_organometallic == False:
                                    centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(
                                    centroid_coords[0], centroid_coords[1], centroid_coords[2]), label='He100{}'.format(str(num_simplified)))
                                    new_component.add_atom(atom)
                                    for comp2 in molecule.components:
                                        if comp2 != component:
                                            new_mol.add_molecule(comp2)
                                    new_mol.add_molecule(new_component)
                                    molecule = new_mol
                                    num_simplified += 1
                                    break
                                    #if 0.85 <= x_inertia/y_inertia <= 1.15:
                                    #    plane = MD.atom_plane(*tuple(a for a in comp.atoms))
                                    #    centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                    #    new_atom_1x = centroid_coords[0] + (0.5 * plane.normal[0])
                                    #    new_atom_1y = centroid_coords[1] + (0.5 * plane.normal[1])
                                    #    new_atom_1z = centroid_coords[2] + (0.5 * plane.normal[2])
                                    #    new_atom_2x = centroid_coords[0] - (0.5 * plane.normal[0])
                                    #    new_atom_2y = centroid_coords[1] - (0.5 * plane.normal[1])
                                    #    new_atom_2z = centroid_coords[2] - (0.5 * plane.normal[2])
                                    #    new_component = ccdc.molecule.Molecule(identifier='centroid')
                                    #    new_atom1 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(new_atom_1x, new_atom_1y, new_atom_1z), label='He100{}'.format(str(num_simplified)))
                                    #    new_atom2 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(new_atom_2x, new_atom_2y, new_atom_2z), label='He101{}'.format(str(num_simplified)))
                                    #    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                    #    na1 = new_component.add_atom(new_atom1)#, new_atom2)
                                    #    na2 = new_component.add_atom(new_atom2)
                                    #    new_component.add_bond(1, na1, na2)
                                    #    for comp2 in molecule.components:
                                    #        if comp2 != component:
                                    #            new_mol.add_molecule(comp2)
                                    #    new_mol.add_molecule(new_component)
                                    #    molecule = new_mol
                                    #    note = note + 'Replaced oblate {} with vector; '.format(comp.formula)
                                    #    break
                                    #else:
                                    #    centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                    #    vector_norm_array = np.array([1,0,0])  #the reorientation moves principle axis to the x-axis
                                    #    rotation_matrix = np.linalg.inv(rotation_matrix)
                                    #    new_vector = rotation_matrix.dot(vector_norm_array)
                                    #    new_atom_1x = centroid_coords[0] + (((box_length-1)/2)  * new_vector[0])
                                    #    new_atom_1y = centroid_coords[1] + (((box_length-1)/2)  * new_vector[1])
                                    #    new_atom_1z = centroid_coords[2] + (((box_length-1)/2)  * new_vector[2])
                                    #    new_atom_2x = centroid_coords[0] - (((box_length-1)/2)  * new_vector[0])
                                    #    new_atom_2y = centroid_coords[1] - (((box_length-1)/2)  * new_vector[1])
                                    #    new_atom_2z = centroid_coords[2] - (((box_length-1)/2)  * new_vector[2])
                                    #    new_component = ccdc.molecule.Molecule(identifier='centroid')
                                    #    new_atom1 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(new_atom_1x, new_atom_1y, new_atom_1z), label='He100{}'.format(str(num_simplified)))
                                    #    new_atom2 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(new_atom_2x, new_atom_2y, new_atom_2z), label='He101{}'.format(str(num_simplified)))
                                    #    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
#
                                    #    na1 = new_component.add_atom(new_atom1)  # , new_atom2)
                                    #    na2 = new_component.add_atom(new_atom2)
                                    #    new_bond = new_component.add_bond(1, na1, na2)
#
                                    #    for comp2 in molecule.components:
                                    #        if comp2 != component:
                                    #            new_mol.add_molecule(comp2)
                                    #    new_mol.add_molecule(new_component)
                                    #    molecule = new_mol
                                    #    note = note + 'Replaced prolate {} with vector; '.format(comp.formula)


        return (molecule, note, num_simplified)

#m_x_contacts is gathered in top.py
def m_x_test(molecule, m_x_contacts, note):
    for cont in m_x_contacts:
        try:
            atom1 = cont[0].label
            atom2 = cont[1].label
            if atom1[-1].isalpha():
                atom1 = atom1[:-1]
            if atom2[-1].isalpha():
                atom2 = atom2[:-1]
            for atom in molecule.atoms:
                if atom.label == atom1:
                    atom1 = atom
                elif atom.label == atom2:
                    atom2 = atom
            molecule.add_bond('Single', atom1, atom2)
        except Exception as e: # should write these out at somepoint to check nothing major is wrong
            pass
            #print(e)
    return(molecule, note)




#------------Pymatgen------------------ functions found online, not sure if they're still be used above

def mol_from_file(fname):
    """
    Reads a file into a pymatgen Molecule. Supported formats include xyz, gaussian,
    and pymatgen JSON. Openbabel is optional but adds additional format options.

    Args:
        fname: the file path string

    Returns:
        a pymatgen Molecule object
    """
    try:
        return Molecule.from_file(fname)
    except:
        print("Error: could not import file " + str(fname) + " to Molecule.\n"
               + "Default supported formats are xyz, gaussian and pymatgen JSON molecules.\n"
               + "Installing openbabel allows for more extensions.")
        return

#def get_box(mol, padding=None):
#    """
#    Given a molecule, find a minimum orthorhombic box containing it.
#    Size is calculated using min and max x, y, and z values,
#    plus the padding defined by the vdw radius
#    For best results, call oriented_molecule first.
#    Args:
#        padding: float (default is 3.4 according to the vdw radius of C)
#    Returns:
#        box: a Box object
#    """
#    #mol, P = reoriented_molecule(self.mol)
#    xyz = mol.cart_coords
#    dims = [0, 0, 0]
#    for i in range(3):
#        dims[i] = np.max(xyz[:, i]) - np.min(xyz[:, i])
#        if padding is not None:
#            dims[i] += padding
#            dims[i] = max([dims[i], 2.0])  # for planar molecules
#        else:
#            ids = np.argsort(xyz[:, i])
#            r = Element(mol[ids[0]].species_string).vdw_radius
#            r += Element(mol[ids[-1]].species_string).vdw_radius
#            dims[i] = max([dims[i] + r, 3.4])  # for planar molecules
#    return Box(dims)

def get_inertia_tensor(mol):
    """
    Calculate the symmetric inertia tensor for a Molecule object. Used to find
    the principal axes of symmetry.

    Args:
        mol: a Molecule object

    Returns:
        a 3x3 numpy array representing a moment of inertia tensor
    """
    mo = mol.get_centered_molecule()
    # Initialize elements of the inertia tensor
    I11 = I22 = I33 = I12 = I13 = I23 = 0.0
    for i in range(len(mo)):
        x, y, z = mo.cart_coords[i]
        m = mo[i].specie.number
        I11 += m * (y ** 2 + z ** 2)
        I22 += m * (x ** 2 + z ** 2)
        I33 += m * (x ** 2 + y ** 2)
        I12 += -m * x * y
        I13 += -m * x * z
        I23 += -m * y * z
    return np.array([[I11, I12, I13],
                  [I12, I22, I23],
                  [I13, I23, I33]])

def reoriented_molecule(mol, nested=False):
    """
    Reorient a molecule so that its principal axes are aligned with the
    identity matrix.

    Args:
        mol: a Molecule object
        nested: internal variable to keep track of how many times the function
            has been called recursively

    Returns:
        new_mol, P: new_mol is a reoriented copy of the original molecule. P is
            the 3x3 rotation matrix used to obtain it.
    """
    def reorient(mol):
        new_mol = mol.get_centered_molecule()
        A = get_inertia_tensor(new_mol)
        #Store the eigenvectors of the inertia tensor
        P = np.transpose(np.linalg.eigh(A)[1])
        if np.linalg.det(P) < 0:
            P[0] *= -1
        #reorient the molecule
        P = SymmOp.from_rotation_and_translation(P,[0,0,0])
        new_mol.apply_operation(P)
        #Our molecule should never be inverted during reorientation.
        if np.linalg.det(P.rotation_matrix) < 0:
            printx("Error: inverted reorientation applied.\n"
            +"(Within reoriented_molecule)", priority=0)
        return new_mol, P
    #If needed, recursively apply reorientation (due to numerical errors)
    iterations = 1
    max_iterations = 100
    new_mol, P = reorient(mol)
    while iterations < max_iterations:
        is_okay = True
        for i in range(3):
            for j in range(3):
                x = np.linalg.eigh(get_inertia_tensor(new_mol))[1][i][j]
                okay = True
                if i == j:
                    #Check that diagonal elements are 0 or 1
                    if (not np.isclose(x, 0)) and (not np.isclose(x, 1)):
                        okay = False
                else:
                    #Check that off-diagonal elements are 0
                    if (not np.isclose(x, 0)):
                        okay = False
        if okay is False:
            #If matrix is not diagonal with 1's and/or 0's, reorient
            new_mol, Q = reorient(new_mol)
            P = Q*P
            iterations += 1
        elif okay is True:
            break
    if iterations == max_iterations:
        printx("Error: Could not reorient molecule after "+str(max_iterations)+" attempts\n"
            +str(new_mol)+"\n"
            +str(get_inertia_tensor(new_mol)), priority=0)
        return False
    return new_mol, P


def get_moment_of_inertia(mol, axis, scale=1.0):
    """
    Calculate the moment of inertia of a molecule about an axis.

    Args:
        mol: a Molecule object
        axis: a 3d axis (list or array) to compute the moment about
        scale: changes the length scale of the molecule

    Returns:
        a scalar value for the moment of inertia about the axis
    """
    #convert axis to unit vector
    axis = axis / np.linalg.norm(axis)
    moment = 0
    for i, a in enumerate(mol):
        v = a.coords
        moment += (scale * np.linalg.norm(np.cross(axis, v)) ) ** 2
    return moment

def test():
    pass

def main():
    remove_water()
