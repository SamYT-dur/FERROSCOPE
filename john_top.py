# script for topological labelling of molecules to find equivalent atoms
# additional code and comments by EDG

# from mercury_interface import MercuryInterface  # NOQA
# 29/4/2020 John change so excluded is a list and can get set by EDG_Master
#
# John started with this EDG_topological_analysis.py from ~15/5/2020 and added extra filters etc
# Needs john_csdtools.py
# 23/5/2020 added exclude C_F_3 type groups
# 27/1/2022 Sam made first draft of simplifying small organics for disorder
# 24/10/2022 Sam added m_x_test for extending metal-halide bond distances so they work with the topology labelling. calls a new function in csdtools

import sys
import math
import os
import glob

from ccdc.io import CrystalReader
from ccdc.io import CrystalWriter
from ccdc.io import EntryReader
from ccdc.molecule import Bond
import john_csdtools
import ccdc.io, ccdc.search, ccdc.entry


excluded = ['H',] #normally set by EDG_Master; included here so EDG_Basic and this script will run standalone
replace_type_by_c = ['U'] #normally set by EDG_Master; included here so EDG_BAsic and this script will run standalone
xy4_exclude = ['ClO4','BrO4','BF4','BH4','IO4'] #normally will be set by EDG_Master; included here so EDG_Basic and this script will run standalone; better to use ABn_list
ABn_list=['C_Cl_3', 'C_F_3', 'C_Br_3'] #normally will be set by EDG_Master; included here so EDG_Basic and this script will run standalone
water_exclude = True
notes_file = 'Outputs/process_notes.txt'
topology_relabel = True

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

    # mercury = MercuryInterface()  # I think this part is unnecessary when running through the terminal, only for Mercury
    # entry_id = mercury.identifier
    #
    # entry = mercury.current_entry
    # crystal = entry.crystal
    # molecule = crystal.molecule
    csd_reader= EntryReader('csd')
    
    entry_id = identifier
    entry = csd_reader.entry(entry_id)
    john_csdtools.entry_id = entry_id
    crystal = entry.crystal

    disordered = crystal.disordered_molecule  # first we take the full molecule, with disorder

    suppressed = [atom for atom in disordered.atoms if atom.label.endswith('?')]  # create a list of suppressed atoms

    #crystal.assign_bonds()

    new_mol = crystal.asymmetric_unit_molecule

    for atom in suppressed:  # here we remove any suppressed atoms from the asymmetric unit
        label = atom.label
        new_mol.remove_atom(new_mol.atom(label))

    molecule = new_mol  # rename new_mol to fit with the rest of the code

#########################################################  M_X_contacts for use in m_x_test. contacts have to be collected at the crystal level and then added in at the molecule level
    halides = ['F', 'Cl', 'Br', 'I', 'O']
    has_metal = False
    has_halide = False

    for atom in molecule.atoms:
        if atom.is_metal:
            has_metal = True
        if atom.atomic_symbol in halides:
            has_halide = True

    if has_metal == True and has_halide == True: # if is_organometallic didn't work for magnesium containing NIWZIE so doing this instead
        contacts = crystal.contacts(intermolecular ='Any', path_length_range=(0, 3), distance_range=(-3.5, 0)) #creates a list of crystallographic contacts. I believe distance_range means contacts up to 3.5 angstroms apart, documentation isn't particularly clear
        #measuring contacts adds a significant amount of processing time.
        metalloids = ['Sb'] # more could be added, e.g. As. only examples in known_ferros contained Sb
        m_x_contacts = []
        for cont in contacts:
            if cont.atoms[0].is_metal or cont.atoms[0].atomic_symbol in metalloids: #next bit checks contact is between halide and metal. doubled up for M-X or X-M bond. probably a cleaner way to do this.
                if cont.atoms[1].atomic_symbol in halides:
                    contact = [cont.atoms[0], cont.atoms[1]]
                    m_x_contacts.append(contact)
            if cont.atoms[1].is_metal or cont.atoms[1].atomic_symbol in metalloids:
                if cont.atoms[0].atomic_symbol in halides:
                    contact = [cont.atoms[0], cont.atoms[1]]
                    m_x_contacts.append(contact)
    else:
        m_x_contacts = []
####################################################################

    # Uncomment this  if you  want bond types assigned
    #molecule = molecule.assign_bond_types(which='unknown')

    waters=0
    tet_anion=0
    note=''


    num_simplified = 0
    #if water_exclude:
    #molecule,waters,note=john_csdtools.remove_water(molecule,waters,note) #take out any water molecules - now done by inertial

    #################### pick and choose what is required for you search here. parameters can be changed inside each function in csd_tools
    molecule=john_csdtools.exclude(molecule,excluded,replace_type_by_c) #exclude and replace atoms


    #molecule, note = john_csdtools.replace_ABn_A(molecule, ABn_list, note) #Sam - i don't think CCl3 removal is working. Td ions seems fine
    #inertial moved before dabco as inertial depends on identifier which can be messed up by dabco
    molecule, note, num_simplified = john_csdtools.simplify_small_organics_inertial(molecule, note,num_simplified)  # Sam's second attempt at simplifaction of small organics using moment of inertia
    molecule, note, num_simplified = john_csdtools.simplify_dabco(molecule, note, num_simplified) #Sam dabco simplification
    #molecule, note = john_csdtools.simplify_small_organics(molecule, note)  # Sam attempt at simplification of small organic molecules - replaced by inertial abov
    molecule, note = john_csdtools.m_x_test(molecule, m_x_contacts, note) #Sam's routine for extending metal-halide bonds
    molecule, note = john_csdtools.sam_Abn_A(molecule, crystal, ABn_list, note)  # sam's replacement for replace_ABn_A #operating on the asym unit does miss a couple entries but not worrying about that for now


    if len(note) > 0:
        print('**Output note for Results2.csv for entry',entry_id,'is:',note)
        with open(notes_file, 'a') as outfile:  # write the note to file
            outfile.write('{0}:{1}\n'.format(entry_id,note))
    if len(note) == 0:
        note='None'
#this bit here does the topology atom relabelling

# Sam - a lot of this stuff still baffles me, I have left it alone as much as possible.
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
    return crystal  # This is all we need when not running through Mercury
    
    # # We must interleave the mol2 and tsv creation because we have a generator!
    # cif_file = os.path.join(mercury.output_directory_path,
    #                          '{}.cif'.format(mercury.output_base))
    # with CrystalWriter(cif_file) as cif_writer:
    #     cif_writer.write(crystal)

        
    #     # Define HTML output file and open it (again we shouldnt need this as it will just print from terminal)
    # html_file = mercury.output_html_file
    # f = open(html_file, "w")
    #
    # tableno = 1
    # # Write report header and title
    # f.write('\n'.join([
    #     '<html>'
    #     '<head>'
    #     '<STYLE TYPE="text/css">',
    #     'p {text-align:justify; font-size:91.67%;}',
    #     'p.ref {text-align:center; font-size:91.67%;}',
    #     'body {font-family:Calibri, Verdana, sans-serif;}',
    #     'table {margin:auto; border-spacing:0px; border-collapse:collapse;}',
    #     'th {padding:2px;}',
    #     'td {text-align:center; padding:2px;}',
    #     'td.ref {text-align:left; padding:2px;}',
    #     '</STYLE>'
    #     '</head>'
    #     '<body>'
    #     '<div>'
    #     '<table>'
    #     '<tr><td class="ref" style="font-size:3em"><b>Crystal Structure Report for %s</b></td>' % entry_id,
    #     '<td><IMG src="%s"></IMG></td></tr>' % mercury.get_ccdc_logo,
    #     '</table>'
    #     '</div>'
    # ]))
    # f.write('<p><center>Wrote file %s</center></p>\n\n' % os.path.normpath(cif_file))

    
def mercury_main():
    """This function does the same thing as above but is used when running from Mercury directly I think.
    None of it is needed"""
    helper = MercuryInterface()
    entry_id = helper.identifier 

    entry = helper.current_entry
    crystal = entry.crystal
    molecule = crystal.molecule
    molecule.assign_bond_types(which='unknown')    


    indexes,strrep = calculate_indexes(molecule)
    
    lookup = []
    for atom in molecule.atoms:
        lookup.append([strrep[atom.index],atom.index])
        
    lookup.sort()    
    
    # From here its just generating an HTML report
    
    # Define HTML output file and open it
    html_file = helper.output_html_file
    f = open(html_file, "w")

    tableno = 1
    # Write report header and title
    f.write('\n'.join([
        '<html>'
        '<head>'
        '<STYLE TYPE="text/css">',
        'p {text-align:justify; font-size:91.67%;}',
        'p.ref {text-align:center; font-size:91.67%;}',
        'body {font-family:Calibri, Verdana, sans-serif;}',
        'table {margin:auto; border-spacing:0px; border-collapse:collapse;}',
        'th {padding:2px;}',
        'td {text-align:center; padding:2px;}',
        'td.ref {text-align:left; padding:2px;}',
        '</STYLE>'
        '</head>'
        '<body>'
        '<div>'
        '<table>'
        '<tr><td class="ref" style="font-size:3em"><b>Crystal Structure Report for %s</b></td>' % entry_id,
        '<td><IMG src="%s"></IMG></td></tr>' % helper.get_ccdc_logo,
        '</table>'
        '</div>'
    ]))
    
    if len(molecule.atoms) == 0:
        f.write('<p><center>This molecule contains no atoms.</center></p>\n\n')
    else:

         
        f.write('\n'.join([
                '<div style="float: middle; text-align: left">'
                '<h2><center>Partitioning</center></h2>']))
        in_table = False
        last_id = ""          
        for ordered in lookup:
            atom = molecule.atoms[ordered[1]]
            if atom.atomic_number == 1: # Skip reporting for H-atoms
                continue
          
            new_table = ordered[0] != last_id

            if new_table and in_table:
                f.write('</table>')

            if new_table:
                f.write('\n'.join([
                '<table align="center"; border="1"; id="cssTable">'
                '<tr><th>Atom</th><th>key</th></tr>'
                ]))
                in_table = True 
                last_id = ordered[0]
                           

            try:
                topo_key = strrep[atom.index]
                f.write('<tr><td>%s</td><td>%s</td></tr>' % (atom.label, topo_key ))
            except TypeError:
                pass
                

    
if __name__ == '__main__':
    #id = input('Enter an entry ID from the CSD:')
#    id = 'FIGNER'
    open(notes_file, 'w').close() #empty the process notes file only needed when running this
    #print('****just deleted notes_file in john_top.py')
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
