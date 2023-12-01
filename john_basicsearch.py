# This code searches through the CSD for each entry subject to some constraints and then analyses any hits with topological analysis
# EDG 28/10/19
# Change number on line 39 to search more/less entries, comment out to search entire database
# Add a hit_list = ['REFCODE'] line to top of top_an() then run that function to get CIFs of specific entries
# 23/4/2020 jsoe added a test for pyroelectric space groups instead of just non-centrosymmetric
# 27/4/2020 jsoe added function to write raw cif files out in case not doing topological analysis
# 29/4/2020 set up so hit_list can get overwritten by master; set up to have variable for highest CSD entry and what to write, added excluded as a list
# 5/5/2022 Sam added write_CSD_cif for use in john_compcif.py
# 3/10/2023 Sam started to clear up and consolidate whole program. Files were backed up before editing. removed hit_list references and let master do all the work

import ccdc.io, ccdc.search, ccdc.entry, ccdc.molecule
import os
import john_top
from ccdc.descriptors import MolecularDescriptors as MD

csd_reader = ccdc.io.EntryReader('CSD')

#pyroelectric point groups
#region
pg_1 = ['1']
pg_1bar = ['2']
pg_2 = ['3', '4', '5']
pg_m = ['6', '7', '8', '9']
pg_2upm = ['10', '11', '12', '13', '14', '15']
pg_222 = ['16', '17', '18', '19', '20', '21', '22', '23', '24']
pg_mm2 = ['25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46']
pg_mmm = ['47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74']
pg_4 = ['75', '76', '77', '78', '79', '80']
pg_4bar = ['81', '82']
pg_4upm = ['83', '84', '85', '86', '87', '88']
pg_422 = ['89', '90', '91', '92', '93', '94', '95', '96', '97', '98']
pg_4mm = ['99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110']
pg_4bar2m = ['111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122']
pg_4upmmm = ['123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142']
pg_3 = ['143', '144', '145', '146']
pg_3bar = ['147', '148']
pg_32 = ['149', '150', '151', '152', '153', '154', '155']
pg_3m = ['156', '157', '158', '159', '160', '161']
pg_3barm = ['162', '163', '164', '165', '166', '167']
pg_6 = ['168', '169', '170', '171', '172', '173']
pg_6bar = ['174']
pg_6upm = ['175', '176']
pg_622 = ['177', '178', '179', '180', '181', '182']
pg_6mm = ['183', '184', '185', '186']
pg_6barm2 = ['187', '188', '189', '190']
pg_6upmmm = ['191', '192', '193', '194']
pg_23 = ['195', '196', '197', '198', '199']
pg_m3bar = ['200', '201', '202', '203', '204', '205', '206']
pg_432 = ['207', '208', '209', '210', '211', '212', '213', '214']
pg_4bar3m = ['215', '216', '217', '218', '219', '220']
pg_m3barm = ['221', '222', '223', '224', '225', '226', '227', '228', '229', '230']
pyropgs = pg_1 + pg_m + pg_2 + pg_mm2 + pg_3 + pg_3m + pg_4 + pg_4mm + pg_6 + pg_6mm
centro_spacegroups = pg_1bar + pg_2upm + pg_mmm + pg_4upm + pg_4upmmm + pg_3bar + pg_3barm + pg_6upm + pg_6upmmm + pg_m3bar + pg_m3barm
non_centro_spacegroups = pg_1 + pg_2 + pg_m + pg_222 + pg_mm2 + pg_4 + pg_4bar + pg_422 + pg_4mm + pg_4bar2m + pg_3 + pg_32 + pg_3m + pg_6 + pg_6bar + pg_622 + pg_6mm + pg_6barm2 + pg_23 + pg_432 + pg_4bar3m
#endregion
#print('Pyro groups are',pyropgs)


excluded = ['H']  # remove the atoms from the cif here before we relabel them for simplicity will get reset by Master
replace_type_by_c = ['U'] #normally set by john_Master; included here so john_BAsic and this script will run standalone
keep_raw_cifs = True
highest_entry = 2000 #normally set by john_master.py
nth_entry = 1 #nomrally set by john_master.py
#john_topological_analysis.excluded = ['H']
attempts = 0
hits = 0
non = 0

if not os.path.exists('TestSearchOutputs/raw_cifs'):  # creates a new dir for all the outputted files
    os.makedirs('TestSearchOutputs/raw_cifs')
out_dir = 'TestSearchOutputs/raw_cifs'

if not os.path.exists('TestSearchOutputs/top_an_cifs'):
    os.makedirs('TestSearchOutputs/top_an_cifs')
top_an_out_dir = 'TestSearchOutputs/top_an_cifs/'
notes_file = 'MasterOutputs/process_notes.txt'

def search():
    """Searches the for non-centrosymmetric pyroelectric crystals and gives an output of all the cifs"""
    attempts = 0
    hits = 0
    non = 0

    search_settings = ccdc.search.Search.Settings()  # defines some search settings, easy to add more/change, see full list in docs
    search_settings.has_3d_coordinates = True
    #search_settings.only_organic = True
    search_settings.no_errors = True
    #search_settings.no_metals = True
    #search_settings.not_polymeric = True
    for i, entry in enumerate(csd_reader):
        attempts = i
        if i % 5000 == 0:  # progress updates
            print('Searched through {} CSD entries of {}.  Number passing symmetry test (probably pyro) found {}'.format(i, len(csd_reader),non))
        if i >= highest_entry:  # for test purposes can just look at first e.g. 2000 entries, set to 1100000 toe search all
            break
        if i % nth_entry == 0: #only write every nth_entry
            if search_settings.test(entry) == True:
                crystal = entry.crystal
                hits += 1
                id = crystal.identifier
                sgold= crystal.spacegroup_symbol
                try:
                    sgnum=str(crystal.spacegroup_number_and_setting[0])
                except: #exception for e.g. EFASCO01,KECYBU15 & MTYHFB03 which doesn't return an sensible space group number
                    sgnum=1
                    print('***Didnt get space group number from csd for',id,'space group',sgold,'setting to 1***')
                if sgnum in pyropgs:
                    pyro = True
                else:
                    pyro= False
                #if crystal.is_centrosymmetric == False:  # we only want the non-centrosymmetric crystals
                if pyro == True:  # we only want the pyroelectric crystals
                    file_name = id + '.cif'
                    output = os.path.join(out_dir, file_name)
                    hit_list.append(id)

                    non += 1
                    #cif = entry.to_string(format='cif')
                    #with open(output, 'w+') as outfile:
                    #    outfile.write(cif)
    print('Excluding preset hit_list, {} search hits from {} attempts. {} are pyroelectric\n'.format(hits, attempts, non))

def run_top_an():
    """This runs the hits from the search through topological analysis for group exclusion and relabelling"""
    open(notes_file, 'w+').close() #empty the process notes file
    #print('**** clearning notes_file')
    cif=0
    for i, hit in enumerate(hit_list):
        cif += 1
        if i % 500 == 0:
            print('Previously processed {} entries from {} through topological analysis'.format(i, len(hit_list)))

        altered = john_top.output_crystal(hit)
        top_an_cif = altered.to_string(format='cif')
        altered_file_name = hit+'_top.cif'
        altered_out = os.path.join(top_an_out_dir, altered_file_name)
        with open(altered_out, 'w+') as outfile:
            outfile.write(top_an_cif)
    print('number of top cifs = '+ str(cif))

def write_CSD_cif():
    """Writes a cif as found in the database to Raw"""
    CSD = ccdc.io.EntryReader('csd')  # now we collect some information on the entry to write to the results
    search_settings = ccdc.search.Search.Settings()  # defines some search settings, easy to add more/change, see full list in docs
    search_settings.has_3d_coordinates = True
    # search_settings.only_organic = True
    search_settings.no_errors = True
    for i, ref_code in enumerate(hit_list):
        entry = CSD.entry(ref_code)
        if search_settings.test(entry) == True:
            crystal = entry.crystal
            id = crystal.identifier
            file_name = id + '.cif'
            output = os.path.join(CSD_out_dir, file_name)

            #disordered = crystal.disordered_molecule  # first we take the full molecule, with disorder

            #suppressed = [atom for atom in disordered.atoms if atom.label.endswith('?')]  # create a list of suppressed atoms
            #new_mol = crystal.asymmetric_unit_molecule

            #for atom in suppressed:  # here we remove any suppressed atoms from the asymmetric unit
            #    label = atom.label
            #    new_mol.remove_atom(new_mol.atom(label))

            #molecule = new_mol  # rename new_mol to fit with the rest of the code
            molecule = crystal.asymmetric_unit_molecule
            for atom in molecule.atoms:
                if atom.label[-1].isalpha():
                    molecule.remove_atom(atom)
                if '_' in atom.label or '?' in atom.label:
                    molecule.remove_atom(atom)
                if atom.atomic_symbol == 'D': #required for reading into pymatgen
                    atom.atomic_symbol = 'H'

                #pairs = MD.AtomDistanceSearch(molecule)
                #print(atom)
                #print(pairs.atoms_within_range(point=atom.coordinates, radius=1))



                #new_atom = ccdc.molecule.Atom(atomic_symbol=atom.atomic_symbol, coordinates=atom.coordinates)


                #print(atom.occupancy)
                #molecule.remove_atom

                #print(atom.molecular_shell(distance_type='actual', distance_range=(0.0,0.001), atom_selection=[atom]))
                #atom.site.occupancy

            crystal.molecule = molecule
            cif = entry.to_string(format='cif')
            with open(output, 'w+') as outfile:
                outfile.write(cif)


def write_raw_cif():
    """Searches the for non-centrosymmetric pyroelectric crystals and gives an output of all the cifs also removes H"""
    CSD = ccdc.io.EntryReader('csd')  # now we collect some information on the entry to write to the results
    search_settings = ccdc.search.Search.Settings()  # defines some search settings, easy to add more/change, see full list in docs
    search_settings.has_3d_coordinates = True
    # search_settings.only_organic = True
    search_settings.no_errors = True
    for i, ref_code in enumerate(hit_list):
        entry = CSD.entry(ref_code)
        if search_settings.test(entry) == True:
            crystal = entry.crystal
            id = crystal.identifier
            file_name = id + '.cif'
            output = os.path.join(out_dir, file_name)

            molecule = crystal.asymmetric_unit_molecule
            for atom in molecule.atoms:
                if atom.atomic_symbol in replace_type_by_c:
                    atom.atomic_symbol = 'C'
                if atom.atomic_symbol in excluded:
                    molecule.remove_atom(atom)
            crystal.molecule = molecule
            cif = entry.to_string(format='cif')
            with open(output, 'w+') as outfile:
                outfile.write(cif)
            if i % 500 == 0:
                print('processed {} entries from {} through topological analysis'.format(i, len(hit_list)))



if __name__ == '__main__':
#    printonscreen()
    search()
    run_top_an()  #calls john_top.py
    write_raw_cif()

