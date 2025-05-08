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

# Sam - remove common solvent molecules from structures # Not sure how reliable!
def solvents(molecule, note, num_simplified):  #
    dcm = Chem.MolFromSmiles('ClCCl')
    tcm = Chem.MolFromSmiles('ClC(Cl)Cl')
    thf = Chem.MolFromSmiles('C1CCOC1')
    ethanol = Chem.MolFromSmiles('CCO')
    methanol = Chem.MolFromSmiles('CO')
    acetonitrile = Chem.MolFromSmiles('CC#N')
    diethyl_ether = Chem.MolFromSmiles('CCOCC')
    acetone = Chem.MolFromSmiles('CC(=O)C')
    dmso = Chem.MolFromSmiles('CS(=S)C')
    toluene = Chem.MolFromSmiles('Cc1ccccc1')
    pyridine = Chem.MolFromSmiles('c1ccncc1')
    water = Chem.MolFromSmiles('C')
    solvents = ['ClCCl', 'ClC(Cl)Cl', 'C1CCCC1', 'CCC', 'CC', 'CC#C', 'CCCCC', 'CC(=C)C', 'CS(=S)C', 'CC1CCCCC1', 'C1CCCCC1', 'C', 'CCC1CCCCC1']


    for component in molecule.components:
        if len(component.atoms) < 9:
            reader = ccdc.io.CrystalReader('CSD')  # the next lines pick out the pure molecule from the database to avoid asymmetric unit only being part of a molecule. Also undoes any simplification.
            reader_molecule = reader.molecule('{}'.format(molecule.identifier))
            for comp in reader_molecule.components:
                #print(comp.atoms)
                if comp.identifier == component.identifier:
                    for atom in comp.atoms:
                        if atom.atomic_symbol in ['N', 'O']:
                            atom.atomic_symbol = 'C'
                        if atom.atomic_symbol in ['H', 'D']:
                            comp.remove_atom(atom)
                    try:
                        smiles = comp.smiles
                        smiles = smiles.replace("[", "")
                        smiles = smiles.replace("]", "")
                        smiles = smiles.replace("c", "C")
                        smiles = smiles.replace("n", "N")
                        smiles = smiles.replace("o", "O")
                        smiles = smiles.replace("s", "S")
                        #print(smiles)
                        if smiles in solvents:
                            #print(smiles)
                            #centroid_coords = MD.atom_centroid(*tuple(a for a in comp.atoms))
                            #new_atomic_symbol = "He"
                            #central_atom_list = MD.AtomDistanceSearch(comp).atoms_within_range(point=centroid_coords,
                            #                                                                   radius=0.5)  # is there an atom within 0.5 of centroid. Is so, should be left as central atom rather than helium.
                            #if central_atom_list:
                            #    central_atom = central_atom_list[0]
                            #    new_atomic_symbol = central_atom.atomic_symbol
                            new_component = ccdc.molecule.Molecule(identifier='centroid')
                            new_mol = ccdc.molecule.Molecule(identifier='{}'.format(molecule.identifier))
                            #atom = ccdc.molecule.Atom(atomic_symbol=new_atomic_symbol, atomic_number=2, coordinates=(
                            #centroid_coords[0], centroid_coords[1], centroid_coords[2]),
                            #                          label='He100{}'.format(str(num_simplified)))
                            #new_component.add_atom(atom)
                            for comp2 in molecule.components:
                                if comp2 != component:
                                    new_mol.add_molecule(comp2)
                            new_mol.add_molecule(new_component)
                            note = note + 'Replaced solvent molecule {} with centroid; '.format(comp.formula)
                            molecule = new_mol
                            num_simplified += 1
                            break
                    except:
                        with open("errors.txt", "a+") as outfile:
                            outfile.write(molecule.identifier)
                            outfile.write("\n")
    return (molecule, note, num_simplified)

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
