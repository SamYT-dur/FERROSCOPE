# Master script that runs all other parts of the CSD search program
# Type "ccdcpython" in shell then "python MASTER.py"

import os
#import sys
import glob
import time
import datetime
import SEARCH
import SYMMETRY_DETECTION
import STRUCTURE
import TOOLS
import sqlite3
import SUPERIMPOSE

start = time.time()
run_version = datetime.datetime.now()
run_version = run_version.strftime("%d" + '/' + "%m" + '/' + "%Y")
print(run_version)

# To clean commandline output of certain errors - #don't think this works as intended
#sys.stderr = object

# Create and define cif directories
# region
# Subroutines may define these differently to be able to use them standalone

with open('errors.txt', 'w+') as outfile: # clean errors file
	outfile.close()

directories = ['OutputCIFs', 'OutputCIFs/Raw', 'OutputCIFs/CSD', 'OutputCIFs/Modified', 'OutputCIFs/Output', 'OutputCIFs/Misc', 'OutputCIFs/Superimpose', 'OutputCIFs/Errant', 'small_mols']
for directory in directories:
	if not os.path.exists(directory):
		os.makedirs(directory)

out_dir = 'OutputCIFs'
raw_out_dir = 'OutputCIFs/Raw'
CSD_out_dir = 'OutputCIFs/CSD'
top_out_dir = 'OutputCIFs/Modified'
other_out_dir = 'OutputCIFs/Misc'
cif_out_dir = 'OutputCIFs/Output'
comp_out_dir = 'OutputCIFs/Superimpose'
errant_out_dir = 'OutputCIFs/Errant'

SYMMETRY_DETECTION.results = 'OutputCIFs/Results.csv'
SYMMETRY_DETECTION.results2 = 'OutputCIFs/Results2.csv'
SYMMETRY_DETECTION.errors = 'OutputCIFs/errors.txt'
SYMMETRY_DETECTION.notes_file = 'OutputCIFs/process_notes.txt'
STRUCTURE.notes_file = 'OutputCIFs/process_notes.txt'
SEARCH.notes_file = 'OutputCIFs/process_notes.txt'
SEARCH.out_dir = raw_out_dir
SEARCH.CSD_out_dir = CSD_out_dir
SEARCH.top_an_out_dir = top_out_dir
SUPERIMPOSE.comp_out_dir = comp_out_dir   #directory used in SUPERIMPOSE

small_mols_directories = ['mol', 'xyz', 'reo', 'spheres', 'not_spheres']
for directory in small_mols_directories:
	if not os.path.exists('small_mols/{}'.format(directory)):
		os.makedirs('small_mols/{}'.format(directory))
	for file in glob.glob('small_mols/{}/*'.format(directory)):
		os.remove(file)

if os.path.isfile('results_database.sqlite') == False:  # Creates results SQLite database file
	open('results_database.sqlite', 'x')
	con = sqlite3.connect('results_database.sqlite')
	cur = con.cursor()
	cur.execute("""CREATE TABLE results (
				ID integer,
				Refcode text,
				LT_SpaceGroup text,
				HT_SpaceGroup text,
				AizuAllowed text,
				Continuous text,
				Irrep text,
				FerroAxes text,
				SearchHistory text,
				Name text,
				Formula text,
				PhaseTransition text,
				LT_SpaceGroup_No integer,
				HT_SpaceGroup_No integer,
				Centrosymmetry text,
				Overall_rss_C real,
				Overall_rss_P real,
				OverallRMS real,
				MaxRMS real,
				MaxType integer,
				MaxOldType text,
				AtomNo integer,
				SymGroupNo integer,
				RecordedTemp text,
				MeltingTemp text,
				Authors text,
				Journal text,
				Volume integer,
				Page integer,
				Year integer,
				DOI text,
				ProcessNotes text,
				Z integer,
				Zprime integer,
				Density text,
				LT_a real,
				LT_b real,
				LT_c real,
				LT_alpha real,
				LT_beta real,
				LT_gamma real,
				LT_volume real,
				HT_a real,
				HT_b real,
				HT_c real,
				HT_alpha real,
				HT_beta real,
				HT_gamma real,
				HT_volume real,
				VolumeRatio real,
				Analogue text,
				Bioactivity text,
				CCDCnumber text,
				Colour text,
				Database text,
				DepositionDate text,
				LTDisorder text,
				Habit text,
				Organic text,
				Polymeric text,
				Powder text,
				Polymorph text,
				Pressure text,
				Rfactor text,
				RadiationSource text,
				Remarks text,
				Solvent text,
				Source text,
				Synonyms text,
				PRIMARY KEY(ID))
				""")
	print('new database created')
else:
	print('database already exists')

# endregion

#### Answer questions from here downwards
# region
# ->1. delete old files set to True or False; delete unchanged cif and inp/out files set to True or False; write out table of shifts set to True or false
delete_old_files = True
delete_unchanged_cif_inp_out = False
output_shifts = True

# ->2. do the CSD search to extract cifs or just run existing files; for preexisting files set up the name for the topolgically labelled or raw cifs here
do_search = True

# ->3. set top_processed to True or False depending on whether you want to run STRUCTURE.py to do topology labelling
# ->3. 23/5/2020 changed things so top_processed means go through STRUCTURE.py
# ->3. STRUCTURE.py itself does various exclusions/group omissions.  relabelling the atoms is just one of the choices.  Set flags for this in #4
top_processed = True  # go through STRUCTURE.py

# ->4. set up atoms/groups to replace by C then exclude here, applied in STRUCTURE and some (currently) to basicsearch if no topological analysis
# ->4. Sam rewrote ABn_a called sam_Abn_a in csdtools due to old script seeing any A connected to any 4B. Simplified FeCl5 for example and catena structures.
replace_type_by_c = ['N', 'O']  # atoms that get changed to C because they'll be similar in size and could disorder
excluded_atoms = ['H', 'D']  # atoms to delete before symmetry analysis etc
ABn_list = ['S_O_4', 'S_C_4', 'P_O_4', 'Cl_C_4', 'P_C_4', 'C_Cl_3', 'C_F_3', 'C_Br_3', 'Br_C_4', 'I_C_4', 'S_C_4', 'Se_C_4', 'B_F_4', 'B_H_4', 'Fe_Cl_4', 'Fe_Cl_5', 'Fe_Br_4', 'Re_C_4', 'Ga_Cl_4', 'Mn_Br_4', 'Mn_Cl_4', 'Zn_Cl_4', 'Cu_Cl_4', 'Co_Cl_4', 'Al_Cl_4', 'Tl_Cl_4', 'Be_Cl_4', 'Cd_Cl_4', 'Zn_Br_4', 'Co_Br_4', 'Al_Br_4', 'Hg_Br_4', 'Ga_Br_4', 'Zn_I_4', 'Hg_I_4', 'In_I_4', 'Al_I_4']  # groups that likely to be disordered so B atoms removed - Sam -expanded by me, there could be more
topology_relabel = True  # this chooses whether to use topology labelling or not
# which structure corrections to apply
STRUCTURE.excluded_check = True	#replace type by C and remove H's
STRUCTURE.sphericity = False #simplify pseudo-spherical molecules # now handled by centroids below
STRUCTURE.MeCyclo = False #simplify methylcyclopentance and mecyclohexane
STRUCTURE.m_x_bonds = False #measure metal-halogen bonds longer than as defined in mercury (up to 3.5 A). Adds significant computational time
STRUCTURE.tds = True #tetrahedrons not simplified by sphericity, handled here

STRUCTURE.solvents = True
TOOLS.sphericity = False
TOOLS.MeCyclo = False
TOOLS.dabco = False
# ->5. set which entries in CSD to search and add to hit_list here; use 2000000 to search everything 0 to just use hit_list defined below.  Can search every nth entry if wanted
SEARCH.highest_entry = 0
SEARCH.nth_entry = 1

# ->6. pre-populate the hit_list here.  If the list is empty, just search database according to section 5 instructions.
# ->6. or can read from a .gcd file of refcodes; or use 32 known FE list; or a 311 test set or enter your own list
with open('REINTSEE.gcd') as infile:
	gcd_list = infile.readlines()
	gcd_list = [x.strip() for x in gcd_list]
SEARCH.hit_list = gcd_list
#SEARCH.hit_list = ['REVNIQ', 'XIFCAR', 'GOFMOD', 'HEZJIG']
print('Preset hit list is', SEARCH.hit_list)

# ->7. process pyroelectric pointgroups or all non-centrosymmetric: edit SEARCH.py around line 117

# ->8. make false to just do the search and skip running findsym
do_run_fs = True

# ->9. set findsym tolerances here; they default otherwise - Sam's tolerance testing confirmed these were the best values
SYMMETRY_DETECTION.lat_tol = 0.25   # latticeTolerance for cell parameters 1.5 optimal
SYMMETRY_DETECTION.atom_tol = 1.0  # atomicPositionTolerance 1.8 optimal
SYMMETRY_DETECTION.atom_rel_tol = 0.999  # atomicPositionMaxTolerance atom tolerance relative to shortest A - A distance if atom_tol is too large i.e. min(atom_tol,atom_rel_tol*shortest_NN_dist) 0.999

# ->10. Set whether to interact with the database or not. Setting to false will still produce an excel file of results - Sam
interact_with_database = False
new_table = False   # if the interact is true and new table is false then the master table will be updated with new hits
new_table_name = 'TestTable2'  # uncomment and name if new_table is true. make sure name doesn't already exist in the database or will crash
#SYMMETRY_DETECTION.new_table_name = new_table_name  # uncomment if new_table is true

# ->11. set whether to generate superimposed cifs of findsym output and input structures. - Sam
compare_cifs = True

# endregion

# Some more file set-up
# region

SEARCH.excluded = excluded_atoms
SEARCH.replace_type_by_c = replace_type_by_c
STRUCTURE.excluded = excluded_atoms
STRUCTURE.replace_type_by_c = replace_type_by_c
STRUCTURE.ABn_list = ABn_list
STRUCTURE.topology_relabel = topology_relabel
SYMMETRY_DETECTION.new_table = new_table
SYMMETRY_DETECTION.interact_with_database = interact_with_database
SYMMETRY_DETECTION.run_version = run_version
SYMMETRY_DETECTION.other_out_dir = other_out_dir
SYMMETRY_DETECTION.cif_out_dir = cif_out_dir
SYMMETRY_DETECTION.errant_out_dir = errant_out_dir
SYMMETRY_DETECTION.delete_unchanged_cif_inp_out = delete_unchanged_cif_inp_out
SYMMETRY_DETECTION.topology_relabel = topology_relabel
SYMMETRY_DETECTION.output_shifts = output_shifts
SYMMETRY_DETECTION.compare_cifs = compare_cifs

if delete_old_files:
	for item in glob.iglob('OutputCIFs/*/*'):
		os.remove(item)

if interact_with_database:
	run_comments = input('\nRun note --->')
	SYMMETRY_DETECTION.run_comments = run_comments

if interact_with_database is True and new_table is True:
	new_table_name = input('\nNew table name ---->')
	SYMMETRY_DETECTION.new_table_name = new_table_name
	con = sqlite3.connect('results_database.sqlite')
	cur = con.cursor()
	cur.execute("""CREATE TABLE {tab} (
				ID integer,
				Refcode text,
				LT_SpaceGroup text,
				HT_SpaceGroup text,
				AizuAllowed text,
				Continuous text,
				Irrep text,
				FerroAxes text,
				SearchHistory text,
				Name text,
				Formula text,
				PhaseTransition text,
				LT_SpaceGroup_No integer,
				HT_SpaceGroup_No integer,
				Centrosymmetry text,
				Overall_rss_C real,
				Overall_rss_P real,
				OverallRMS real,
				MaxRMS real,
				MaxType integer,
				MaxOldType text,
				AtomNo integer,
				SymGroupNo integer,
				RecordedTemp text,
				MeltingTemp text,
				Authors text,
				Journal text,
				Volume integer,
				Page integer,
				Year integer,
				DOI text,
				ProcessNotes text,
				Z integer,
				Zprime integer,
				Density text,
				LT_a real,
				LT_b real,
				LT_c real,
				LT_alpha real,
				LT_beta real,
				LT_gamma real,
				LT_volume real,
				HT_a real,
				HT_b real,
				HT_c real,
				HT_alpha real,
				HT_beta real,
				HT_gamma real,
				HT_volume real,
				VolumeRatio real,
				Analogue text,
				Bioactivity text,
				CCDCnumber text,
				Colour text,
				Database text,
				DepositionDate text,
				LTDisorder text,
				Habit text,
				Organic text,
				Polymeric text,
				Powder text,
				Polymorph text,
				Pressure text,
				Rfactor text,
				RadiationSource text,
				Remarks text,
				Solvent text,
				Source text,
				Synonyms text,
				PRIMARY KEY(ID))
				""".format(tab = new_table_name))
	con.commit()
	con.close()

# endregion

# Start the search

SEARCH.search()

# Run topology relabelling and structure modifications

if top_processed:
	print("Running cifs through topological analysis")
	if do_search:
		SEARCH.run_top_an()
		SEARCH.write_CSD_cif()  # writes out cifs as found in the database for use in SUPERIMPOSE.py - sam
	SYMMETRY_DETECTION.cif_dir = sorted(glob.glob(top_out_dir+'/*.cif'))  # redefining some more directories
	SYMMETRY_DETECTION.top_processed = True
else:
	print("Will run raw cifs without topological analysis")
	if do_search:
		SEARCH.write_raw_cif() # option to write the raw cifs i.e. before topological analysis
	SYMMETRY_DETECTION.cif_dir = sorted(glob.glob(raw_out_dir+'/*.cif'))  # redefining some more directories
	SYMMETRY_DETECTION.top_processed = False

# Run the structures through the symmetry detection

if do_run_fs:
	SYMMETRY_DETECTION.main()

end = time.time()

#with open('OutputCIFS/Results.csv', 'a') as outfile:
	#outfile.write('Script took '+str(end-start)+' seconds\n')

print('Script took '+str(end-start)+' seconds')


