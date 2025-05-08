import math
import os
from iodata import load_one, dump_one
import ccdc.molecule
import ccdc.crystal
import ccdc.search
import ccdc.io
#import ccdc.conformer
from ccdc.descriptors import MolecularDescriptors as MD
from ccdc.descriptors import GeometricDescriptors as GD

import traceback
import sys
import numpy as np
import sqlite3
# Imports
# ------------------------------
# External Libraries
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.symmetry.analyzer import generate_full_symmops

# PyXtal imports
from pyxtal.operations import *
from pyxtal.molecule import *
from pyxtal.viz import *
from pyxtal.database.collection import Collection

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, rdmolfiles
from rdkit import DataStructs

# imports for globularity testing
from scipy.spatial import ConvexHull
import miniball
from math import pi, sqrt

# File set-up for running standalone. Otherwise, these are set by MASTER.py
# region
replace_type_by_c = ['N', 'O']  # atoms that get changed to C because they'll be similar in size and could disorder
excluded_atoms = ['H', 'D']  # atoms to delete before symmetry analysis etc
ABn_list = ['C_Cl_3', 'C_F_3', 'C_Br_3', 'Cl_O_4', 'Br_O_4', 'I_O_4', 'S_O_4', 'Cl_C_4', 'Br_C_4', 'I_C_4', 'B_F_4', 'B_H_4', 'Fe_Cl_4', 'Fe_Cl_5', 'Fe_Br_4', 'Re_O_4', 'Re_C_4', 'Ga_Cl_4', 'Mn_Br_4', 'Mn_Cl_4', 'Zn_Cl_4', 'Cu_Cl_4', 'Co_Cl_4', 'Al_Cl_4', 'Tl_Cl_4', 'Be_Cl_4', 'Cd_Cl_4', 'Zn_Br_4', 'Co_Br_4', 'Al_Br_4', 'Hg_Br_4', 'Ga_Br_4', 'Zn_I_4', 'Hg_I_4', 'In_I_4', 'Al_I_4']  # groups that likely to be disordered so B atoms removed - Sam -expanded by me, there could be more
notes_file = 'Outputs/process_notes.txt'
# endregion


def exclude(molecule, excluded, replace_type_by_c):
    for atom in molecule.atoms:
        if atom.atomic_symbol in replace_type_by_c:
            atom.atomic_symbol = 'C'
        if atom.atomic_symbol in excluded:
            molecule.remove_atom(atom)
    return molecule


# Sam's attempt at replacing replace_ABn_a. fixes a lot of issues with the previous iteration but has problems of its own. Operating on the asymmetric unit means that some tetrahedrons will look like A-B2 rather than A-B4 for example and won't be simplified.
def sam_Abn_A(molecule, crystal, ABn_list, note):
    to_remove = []
    crystmol = crystal.molecule
    crystmol.remove_hydrogens()
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


# Sam's method for simplifying 18-crown-6. Only gains two entries and loses POXQUR as one crown is on a  special position. Currently not used in settings
def simplify_crown(molecule, note, num_simplified):
    crown = Chem.MolFromSmiles('C1CCCCCCCCCCCCCCCCC1')
    for component in molecule.components:
        print(component.atoms)
        if len(component.atoms) < 20:
            smiles = component.smiles
            print(smiles)
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                overlap = mol.GetSubstructMatch(crown)
                if overlap != ():
                    print('yes')
                    plane = MD.atom_plane(*tuple(a for a in component.atoms))
                    centroid_coords = MD.atom_centroid(*tuple(a for a in component.atoms))
                    new_atom_1x = centroid_coords[0] + (0.5 * plane.normal[0])
                    new_atom_1y = centroid_coords[1] + (0.5 * plane.normal[1])
                    new_atom_1z = centroid_coords[2] + (0.5 * plane.normal[2])
                    new_atom_2x = centroid_coords[0] - (0.5 * plane.normal[0])
                    new_atom_2y = centroid_coords[1] - (0.5 * plane.normal[1])
                    new_atom_2z = centroid_coords[2] - (0.5 * plane.normal[2])
                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                    new_atom1 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2,
                                                   coordinates=(new_atom_1x, new_atom_1y, new_atom_1z),
                                                   label='He100{}'.format(str(num_simplified)))
                    new_atom2 = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2,
                                                   coordinates=(new_atom_2x, new_atom_2y, new_atom_2z),
                                                   label='He101{}'.format(str(num_simplified)))
                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                    na1 = new_component.add_atom(new_atom1)  # , new_atom2)
                    na2 = new_component.add_atom(new_atom2)
                    new_component.add_bond(1, na1, na2)
                    for comp2 in molecule.components:
                        if comp2 != component:
                            new_mol.add_molecule(comp2)
                    new_mol.add_molecule(new_component)
                    molecule = new_mol
                    num_simplified += 1
                    note = note + 'simplified crown'
                    #break
                else:
                    pass
                    # print('no')
            except Exception as e:
                print(e)
                pass
    return molecule, note, num_simplified


# Sam - straightforward removal of any molecule containing a DABCO-like substructure.
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
                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2,
                                              coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]),
                                              label='He100{}'.format(str(num_simplified)))
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


# Sam - removes items with methyl-cyclopentane/methyl-cyclohexane substructure
def simplify_MeCyclo(molecule, note, num_simplified):
    MeC5 = Chem.MolFromSmiles('CC1CCCC1')
    MeC6 = Chem.MolFromSmiles('CC(CCC1)CC1')
    for component in molecule.components:
        if len(component.atoms) < 9:
            smiles = component.smiles
            try:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                overlap = mol.GetSubstructMatch(MeC5)
                overlap2 = mol.GetSubstructMatch(MeC6)
                if overlap != () or overlap2 != ():
                    centroid_coords = MD.atom_centroid(*tuple(a for a in component.atoms))
                    new_component = ccdc.molecule.Molecule(identifier='centroid')
                    new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                    atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2,
                                              coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]),
                                              label='He100{}'.format(str(num_simplified)))
                    new_component.add_atom(atom)
                    for comp2 in molecule.components:
                        if comp2 != component:
                            new_mol.add_molecule(comp2)
                    new_mol.add_molecule(new_component)
                    molecule = new_mol
                    num_simplified += 1
            except TypeError: # fails if smiles cannot be made
                continue
            except Exception as e:
                with open('errors.txt', 'a+') as outfile:
                    outfile.write(str(molecule.identifier) + ' : TOOLS.simplify_MeCyclo: ' + str(e) + '\n')
                continue
    return molecule, note, num_simplified

# Sam - moment of inertia based organic simplification
def centroids(molecule, note, num_simplified):  # Sam's V2
    MeC5 = Chem.MolFromSmiles('CC1CCCC1')
    MeC6 = Chem.MolFromSmiles('CC(CCC1)CC1')
    dabco = Chem.MolFromSmiles('C1CC2CCC1CC2')

    for component in molecule.components:
        if len(component.atoms) < 14:
            reader = ccdc.io.CrystalReader('CSD')  # the next lines pick out the pure molecule from the database to avoid asymmetric unit only being part of a molecule. Also undoes any simplification to enable more accurate inertia values.
            reader_molecule = reader.molecule('{}'.format(molecule.identifier))
            carbons = 0
            other_atoms = 0
            for atom in component.atoms:
                atom_remove = ['F', 'Cl', 'Br', 'I', 'P', 'O', 'N']
                if atom.atomic_symbol == 'C':
                    carbons += 1
                if atom.atomic_symbol in atom_remove:
                    other_atoms += 1
            if carbons + other_atoms == len(component.atoms) and other_atoms != len(component.atoms):
                for comp in reader_molecule.components:
                    if comp.identifier == component.identifier:
                        try:
                            for atom in comp.atoms: #remove dueterium atoms
                                label = atom.label
                                if label.startswith('D') and label[1].isdigit():
                                    comp.remove_atom(atom)
                            comp.add_hydrogens(mode='missing', add_sites=True)
                        except RuntimeError:
                            for atom in comp.atoms:  # prevents 'atom has no coordinates error'
                                if atom.coordinates == None:
                                    comp.remove_atom(atom)
                            pass
                        except Exception as e:
                            with open('errors.txt', 'a+') as outfile:
                                outfile.write(
                                    str(molecule.identifier) + ' : TOOLS.sphericity(add hydrogens): ' + str(e) + '\n')
                        mol_file = comp.to_string(
                            format='mol2')  # csd doesn't output .xyz files so output and write out mol2 file
                        if os.path.isfile('small_mols/mol/{}{}.mol2'.format(reader_molecule.identifier,
                                                                            str(comp.identifier))) is False:
                            with open('small_mols/mol/{}{}.mol2'.format(reader_molecule.identifier,
                                                                        str(comp.identifier)), 'w+') as datafile:
                                datafile.write(mol_file)
                        else:
                            with open('small_mols/mol/{}{}.mol2'.format(reader_molecule.identifier,
                                                                        str(comp.identifier)), 'w+') as datafile:
                                datafile.write(mol_file)
                        try:
                            iodata_mol = load_one('small_mols/mol/{}{}.mol2'.format(reader_molecule.identifier,
                                                                                    str(comp.identifier)))  # these two lines us iodata module to convert mol2 to xyz and write out
                            dump_one(iodata_mol, 'small_mols/xyz/{}{}.xyz'.format(reader_molecule.identifier,
                                                                                  str(comp.identifier)))
                            pymatgen_mol = mol_from_file(
                                'small_mols/xyz/{}{}.xyz'.format(reader_molecule.identifier,
                                                                 str(comp.identifier)))
                            reo_mol = reoriented_molecule(pymatgen_mol)[0]
                            rotation_matrix = reoriented_molecule(pymatgen_mol)[1].rotation_matrix
                            reo_mol.to(filename='small_mols/reo/{}.json'.format(reader_molecule.identifier))
                            x_inertia = get_moment_of_inertia(reo_mol, [1, 0, 0])
                            y_inertia = get_moment_of_inertia(reo_mol, [0, 1, 0])
                            z_inertia = get_moment_of_inertia(reo_mol, [0, 0, 1])
                            spherocity_index = (3 * x_inertia) / (x_inertia + y_inertia + z_inertia)  # with weights

                            with open('inertias.txt', 'a+') as outfile:
                                outfile.write(molecule.identifier)
                                outfile.write('\n')
                                outfile.write('{} {} {}'.format(x_inertia, y_inertia, z_inertia))
                            posns = []

                            for atom in reo_mol.sites:
                                atom_coord = []
                                for coord in atom.coords:
                                    atom_coord.append(coord)
                                posns.append(atom_coord)
                        except TypeError: # can't reorient molecule leading to False molecule
                            break
                        except Exception as e:
                            with open('errors.txt', 'a+') as outfile:
                                outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(inertia_calc): ' + str(e) + '\n')
                            break
                        if x_inertia < 3 and y_inertia < 3 and z_inertia < 5:  # removing water and ammonium
                            new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                            for comp2 in molecule.components:
                                if comp2 != component:
                                    new_mol.add_molecule(comp2)
                            molecule = new_mol
                            note = note + 'Removed {}; '.format(comp.formula)
                            break

                        for atom in comp.atoms:
                            if atom.atomic_symbol in ['N', 'O']:
                                atom.atomic_symbol = 'C'
                            if atom.atomic_symbol in ['H', 'D']:
                                comp.remove_atom(atom)
                        smiles = comp.smiles
                        overlap = ()
                        overlap2 = ()
                        overlap3 = ()
                        try:
                            mol = Chem.MolFromSmiles(smiles, sanitize=False)
                            overlap = mol.GetSubstructMatch(dabco)
                            if len(comp.atoms) < 9:
                                overlap2 = mol.GetSubstructMatch(MeC5)
                                overlap3 = mol.GetSubstructMatch(MeC6)
                        except:
                            pass


                        if (0.9 <= spherocity_index <= 1.1) or overlap != () or overlap2 != () or overlap3 != ():
                            centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                            new_atomic_symbol = "He"
                            central_atom_list = MD.AtomDistanceSearch(comp).atoms_within_range(point=centroid_coords,
                                                                                               radius=0.5)  # is there an atom within 0.5 of centroid. Is so, should be left as central atom rather than helium.
                            if central_atom_list:
                                central_atom = central_atom_list[0]
                                new_atomic_symbol = central_atom.atomic_symbol

                            new_component = ccdc.molecule.Molecule(identifier='centroid')
                            new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                            atom = ccdc.molecule.Atom(atomic_symbol=new_atomic_symbol, atomic_number=2, coordinates=(
                            centroid_coords[0], centroid_coords[1], centroid_coords[2]),
                                                      label='He100{}'.format(str(num_simplified)))
                            new_component.add_atom(atom)
                            for comp2 in molecule.components:
                                if comp2 != component:
                                    new_mol.add_molecule(comp2)
                            new_mol.add_molecule(new_component)
                            note = note + 'Replaced spherical molecule {} with centroid; '.format(comp.formula)
                            molecule = new_mol
                            num_simplified += 1
                            try:
                                reo_mol.to(filename='small_mols/spheres/{}.xyz'.format(small_organics_molecule.identifier),
                                       fmt='xyz')
                            except:
                                pass
                            break
    return (molecule, note, num_simplified)

# Sam - moment of inertia based organic simplification
def simplify_small_organics_inertial(molecule, note, num_simplified):  # Sam's V2
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
            if carbons + other_atoms == len(component.atoms) and other_atoms != len(component.atoms):
                small_organics_reader = ccdc.io.CrystalReader('CSD')  # the next lines pick out the pure molecule from the database to avoid asymmetric unit only being part of a molecule. Also undoes any simplification to enable more accurate inertia values.
                small_organics_molecule = small_organics_reader.molecule('{}'.format(molecule.identifier))

                for comp in small_organics_molecule.components:
                    if comp.identifier == component.identifier:
                        try:  # incase adding hydrogens fails e.g. unknown bond type
                            comp.add_hydrogens(mode='missing', add_sites=True)
                        except RuntimeError:
                            for atom in comp.atoms:  # prevents 'atom has no coordinates error'
                                if atom.coordinates == None:
                                    comp.remove_atom(atom)
                            pass
                        except Exception as e:
                            with open('errors.txt', 'a+') as outfile:
                                outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(add hydrogens): ' + str(e) + '\n')
                        mol_file = comp.to_string(format='mol2')  # csd doesn't output .xyz files so output and write out mol2 file
                        if os.path.isfile('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier))) is False:
                            with open('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)), 'w+') as datafile:
                                datafile.write(mol_file)
                        else:
                            with open('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)), 'w+') as datafile:
                                datafile.write(mol_file)
                        try:
                            iodata_mol = load_one('small_mols/mol/{}{}.mol2'.format(small_organics_molecule.identifier, str(comp.identifier)))  # these two lines us iodata module to convert mol2 to xyz and write out
                            dump_one(iodata_mol, 'small_mols/xyz/{}{}.xyz'.format(small_organics_molecule.identifier, str(comp.identifier)))
                            pymatgen_mol = mol_from_file('small_mols/xyz/{}{}.xyz'.format(small_organics_molecule.identifier,str(comp.identifier)))
                            reo_mol = reoriented_molecule(pymatgen_mol)[0]
                            rotation_matrix = reoriented_molecule(pymatgen_mol)[1].rotation_matrix
                            reo_mol.to(filename='small_mols/reo/{}.json'.format(small_organics_molecule.identifier))
                            x_inertia = get_moment_of_inertia(reo_mol, [1, 0, 0])
                            y_inertia = get_moment_of_inertia(reo_mol, [0, 1, 0])
                            z_inertia = get_moment_of_inertia(reo_mol, [0, 0, 1])
                            spherocity_index = (3 * x_inertia) / (x_inertia + y_inertia + z_inertia)  # with weights

                            with open('inertias.txt', 'a+') as outfile:
                                outfile.write(molecule.identifier)
                                outfile.write('\n')
                                outfile.write('{} {} {}'.format(x_inertia,y_inertia,z_inertia))
                            posns = []

                            for atom in reo_mol.sites:
                                atom_coord = []
                                for coord in atom.coords:
                                    atom_coord.append(coord)
                                posns.append(atom_coord)
                            posns2 = np.array(posns)
                            if (len(component.atoms)) > 3:
                                try:
                                    globularity = sphericity(posns2, core=False)
                                except Exception as error:
                                    globularity = 1  # if globularity fails, just call it spherical and only sphericity will be accounted for
                                    if "QH6154 Qhull precision error" in str(error):
                                        continue
                                    else:
                                        with open('errors.txt', 'a+') as outfile:
                                            outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(globularity): ' + str(error) + '\n')
                            else:
                                globularity = 1

                        except TypeError: # can't reorient molecule leading to False molecule
                            break
                        except Exception as e:
                            with open('errors.txt', 'a+') as outfile:
                                outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(inertia_calc): ' + str(e) + '\n')
                            break

                        if x_inertia < 3 and y_inertia < 3 and z_inertia < 5:  # removing water and ammonium
                            new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                            for comp2 in molecule.components:
                                if comp2 != component:
                                    new_mol.add_molecule(comp2)
                            molecule = new_mol
                            note = note + 'Removed {}; '.format(comp.formula)
                            break

                        if (0.78 <= spherocity_index <= 1.1) and globularity > 0.05:  # simplify spherical molecules by replacing with a helium at the centroid. globularity check is to prevent SOMDOO being mistaken
                            try:
                                centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                                new_component = ccdc.molecule.Molecule(identifier='centroid')
                                new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                                atom = ccdc.molecule.Atom(atomic_symbol='He', atomic_number=2, coordinates=(centroid_coords[0], centroid_coords[1], centroid_coords[2]),label='He100{}'.format(str(num_simplified)))
                                new_component.add_atom(atom)
                                for comp2 in molecule.components:
                                    if comp2 != component:
                                        new_mol.add_molecule(comp2)
                                new_mol.add_molecule(new_component)
                                note = note + 'Replaced spherical molecule {} with point; '.format(comp.formula)
                                molecule = new_mol
                                num_simplified += 1
                                reo_mol.to(filename='small_mols/spheres/{}.xyz'.format(small_organics_molecule.identifier), fmt='xyz')
                                break
                            except Exception as e:
                                with open('errors.txt', 'a+') as outfile:
                                    outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(centroid): ' + str(e) + '\n')
                        else:
                            reo_mol.to(filename='small_mols/not_spheres/{}.xyz'.format(small_organics_molecule.identifier),
                                       fmt='xyz')

    return (molecule, note, num_simplified)


# m_x_contacts is gathered in top.py
def m_x_test(molecule, crystal, note):
    halides = ['F', 'Cl', 'Br', 'I', 'O']
    has_metal = False
    has_halide = False

    for atom in molecule.atoms:
        if atom.is_metal:
            has_metal = True
        if atom.atomic_symbol in halides:
            has_halide = True

    if has_metal == True and has_halide == True:  # if is_organometallic didn't work for magnesium containing NIWZIE so doing this instead
        contacts = []
        try:
            contacts = crystal.contacts(intermolecular='Any', path_length_range=(0, 3), distance_range=(-3.5, 0))  # creates a list of crystallographic contacts. I believe distance_range means contacts up to 3.5 angstroms apart, documentation isn't particularly clear
        except RuntimeError: # bond already exists
            pass
        except Exception as e:
            with open('errors.txt', 'a+') as outfile:
                outfile.write(str(molecule.identifier) + ' : TOOLS.m_x(read): ' + str(e) + '\n')
            pass
        # measuring contacts adds a significant amount of processing time.
        metalloids = ['Sb']  # more could be added, e.g. As. only examples in known_ferros contained Sb
        m_x_contacts = []
        for cont in contacts:
            if cont.atoms[0].is_metal or cont.atoms[0].atomic_symbol in metalloids:  # next bit checks contact is between halide and metal. doubled up for M-X or X-M bond. probably a cleaner way to do this.
                if cont.atoms[1].atomic_symbol in halides:
                    contact = [cont.atoms[0], cont.atoms[1]]
                    m_x_contacts.append(contact)
            if cont.atoms[1].is_metal or cont.atoms[1].atomic_symbol in metalloids:
                if cont.atoms[0].atomic_symbol in halides:
                    contact = [cont.atoms[0], cont.atoms[1]]
                    m_x_contacts.append(contact)
    else:
        m_x_contacts = []

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
        except RuntimeError: # bond already exists
            continue
        except Exception as e:
            with open('errors.txt', 'a+') as outfile:
                outfile.write(str(molecule.identifier) + ' : TOOLS.m_x(write): ' + str(e) + '\n')
            pass
            # print(e)KIGSAV
    return (molecule, note)


# ------------Pymatgen------------------ functions found online, leaving here in case used above

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
    except Exception as e:
        #print("Error: could not import file " + str(fname) + " to Molecule.\n"
        #      + "Default supported formats are xyz, gaussian and pymatgen JSON molecules.\n"
        #      + "Installing openbabel allows for more extensions.")
        with open('errors.txt', 'a+') as outfile:
            outfile.write(str(molecule.identifier) + ' : TOOLS.sphericity(molfromfile): ' + str(e) + '\n')
        return

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
        # Store the eigenvectors of the inertia tensor
        P = np.transpose(np.linalg.eigh(A)[1])
        if np.linalg.det(P) < 0:
            P[0] *= -1
        # reorient the molecule
        P = SymmOp.from_rotation_and_translation(P, [0, 0, 0])
        new_mol.apply_operation(P)
        # Our molecule should never be inverted during reorientation.
        if np.linalg.det(P.rotation_matrix) < 0:
            printx("Error: inverted reorientation applied.\n"
                   + "(Within reoriented_molecule)", priority=0)
        return new_mol, P

    # If needed, recursively apply reorientation (due to numerical errors)
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
                    # Check that diagonal elements are 0 or 1
                    if (not np.isclose(x, 0)) and (not np.isclose(x, 1)):
                        okay = False
                else:
                    # Check that off-diagonal elements are 0
                    if (not np.isclose(x, 0)):
                        okay = False
        if okay is False:
            # If matrix is not diagonal with 1's and/or 0's, reorient
            new_mol, Q = reorient(new_mol)
            P = Q * P
            iterations += 1
        elif okay is True:
            break
    if iterations == max_iterations:
        printx("Error: Could not reorient molecule after " + str(max_iterations) + " attempts\n"
               + str(new_mol) + "\n"
               + str(get_inertia_tensor(new_mol)), priority=0)
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
    # convert axis to unit vector
    axis = axis / np.linalg.norm(axis)
    moment = 0
    for i, a in enumerate(mol):
        v = a.coords
        moment += (scale * np.linalg.norm(np.cross(axis, v))) ** 2
    return moment


# -------------------------- Globularity function
def sphericity(molecule, core=True, return_components=False):
    """
    Calculate the sphericity of a molecule.

    Parameters:
    -----------
    molecule: object
        The molecule object.
    core: bool, optional (default=True)
        If True, only the core part of the molecule will be considered for sphericity calculation.
    return_components: bool, optional (default=False)
        If True, returns the components of the sphericity calculation: volumes of the convex hull
        and the bounding sphere.

    Returns:
    --------
    float or tuple:
        The sphericity value, or a tuple containing the volumes of the convex hull and the bounding
        sphere if return_components is True.
    """

    # if not mt.is_connected(molecule, wrapping=False):
    #    raise ValueError("Error: molecule is not connected without wrapping")

    def calculate_sphericity(coords):
        """
        Helper function to calculate sphericity from the input coordinates.

        Parameters:
        -----------
        coords: array_like
            The coordinates of the atoms in the molecule.

        Returns:
        --------
        tuple:
            The volumes of the convex hull and the bounding sphere.
        """
        CH = ConvexHull(coords)
        V_conv = CH.volume
        C, r2 = miniball.get_bounding_ball(coords)
        V_ball = 4 * pi * sqrt(r2) ** 3 / 3
        return V_conv, V_ball

    if core:
        adj_mat = mt.adjacency_matrix(molecule, cutoff=1.1, diag=0)
        G = nx.from_numpy_matrix(adj_mat)
        bridge = nx.bridges(G)
        G.remove_edges_from(bridge)
        G.remove_nodes_from(list(nx.isolates(G)))

        if len(list(G)) > 3:
            atoms = molecule[list(G)]
            mols = mt.decompose(atoms)

            V_convs = []
            V_balls = []
            for mol in mols:
                Rs = mol.positions
                V_conv, V_ball = calculate_sphericity(Rs)
                V_convs.append(V_conv)
                V_balls.append(V_ball)

            if return_components:
                return V_convs, V_balls

            return max(V_conv / V_ball for V_conv, V_ball in zip(V_convs, V_balls))

    else:
        Rs = molecule#.positions
        V_conv, V_ball = calculate_sphericity(Rs)
        if return_components:
            return V_conv, V_ball

        return V_conv / V_ball
