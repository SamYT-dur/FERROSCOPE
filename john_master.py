# Master script that runs all other parts of the CSD search program
# Type "ccdcpython" in shell then "python john_master.py"
#
# JSOE notes
# 27/4/2020 added option of true/false for topological analysis
# 29/4/2020 set things up so could control most things from master without having to change other scripts
# 29/4/2020 i.e. you just have to answer questions ->1 to ->9 below and there are minimal edits in other .py files
# 22/5/2020 set up more possibilities to exclude water etc; uses john_top.py and john_csdtools.py
# 27/5/2020 once extra searches and shifts included saved as john_master.py
# 27/5/2020 will call john_basicsearch, john_top and john_runfs, some of these need john_csdtools.py to work
# 28/1/2022 Sam updated findsym to 7.1.3
# 25/2/2022 Sam couldnt figure out why some top.cifs were getting deleted (something in runfs) so just wrote out a second folder
# 28/4/2022 Sam added compcif prototype. currently running through external john_compcif.py
# 5/5/2022  Sam added new folder to MasterOutputs called CSD cifs for writing out cifs as seen in the database for use in john_compcif.py
# 3/10/2023 Sam started to clear up and consolidate whole program. Files were backed up before editing. deleted references to inertial_database here and in csd_tools


import os
import glob
import time
import datetime
import john_basicsearch
import john_runfs
import john_top
import sqlite3
import john_compcif3

start = time.time()
run_version = datetime.datetime.now()
run_version = run_version.strftime("%d" + '/' + "%m" + '/' + "%Y")
print(run_version)

# set up directories etc here.  Keep the same basic names to avoid confusion (inherited from Elliot in part)
# definitions here will override definitions in the other .py files; done like this so individual .py's will run stand-alone
if not os.path.exists('MasterOutputs'):  # creates a new dir for all the outputted files
	os.makedirs('MasterOutputs')
out_dir = 'MasterOutputs'

if not os.path.exists('MasterOutputs/Raw_cifs'):  # creates a new dir for all the raw cifs
	os.makedirs('MasterOutputs/Raw_cifs')
raw_out_dir = 'MasterOutputs/Raw_cifs'

if not os.path.exists('MasterOutputs/CSD_cifs'):  # creates a new dir for all the raw cifs - sam
	os.makedirs('MasterOutputs/CSD_cifs')
CSD_out_dir = 'MasterOutputs/CSD_cifs'

if not os.path.exists('MasterOutputs/Top'):  # creates a new dir for all the top_an cifs
	os.makedirs('MasterOutputs/Top')
top_out_dir = 'MasterOutputs/Top'

if not os.path.exists('MasterOutputs/Other'):  # creates an output dir for all the inp/out files
	os.makedirs('MasterOutputs/Other')
other_out_dir = 'MasterOutputs/Other'

if not os.path.exists('MasterOutputs/Findsym_cifs'):  # a new dir for just the outputted findsym_cif files
	os.makedirs('MasterOutputs/Findsym_cifs')
cif_out_dir = 'MasterOutputs/Findsym_cifs'

if not os.path.exists('MasterOutputs/Compare_cifs'):  # Sam - a new dir for just the outputted findsym_cif files
	os.makedirs('MasterOutputs/Compare_cifs')
comp_out_dir = 'MasterOutputs/Compare_cifs'


if not os.path.exists('MasterOutputs/Errant'):  # a new dir for just the errant_cif files for tracking
	os.makedirs('MasterOutputs/Errant')

files = glob.glob('small_mols/mol/*') #Sam - this chunk just empties the diectories for small_organics
files2 = glob.glob('small_mols/xyz/*')
for f in files:
	os.remove(f)
for f in files2:
	os.remove(f)

if os.path.isfile('/home/dpqk92/PycharmProjects/csd_search_john/results_database.sqlite') == False:  #sam - creates database file with one master table called results the first time the program is run
	open('/home/dpqk92/PycharmProjects/csd_search_john/results_database.sqlite', 'x')
	con = sqlite3.connect('/home/dpqk92/PycharmProjects/csd_search_john/results_database.sqlite')
	cur = con.cursor()
	cur.execute("""CREATE TABLE results (
				ID integer,
				Refcode text,
				LT_SpaceGroup text,
				HT_SpaceGroup text,
				AizuAllowed text,
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

errant_out_dir = 'MasterOutputs/Errant'
john_runfs.results = 'MasterOutputs/Results.csv'
john_runfs.results2 = 'MasterOutputs/Results2.csv'
john_runfs.errors = 'MasterOutputs/errors.txt'
john_runfs.notes_file = 'MasterOutputs/process_notes.txt'
john_top.notes_file = 'MasterOutputs/process_notes.txt'
john_basicsearch.notes_file = 'MasterOutputs/process_notes.txt'
john_basicsearch.out_dir = raw_out_dir
john_basicsearch.CSD_out_dir = CSD_out_dir
john_basicsearch.top_an_out_dir = top_out_dir
john_compcif3.comp_out_dir = comp_out_dir # Sam - directory used in compcif3


#### Answer questions from here downwards

# ->1. delete old files set to True or False; delete unchanged cif and inp/out files set to True or False; write out table of shifts set to True or false
delete_old_files = True
delete_unchanged_cif_inp_out = True
output_shifts = True

# ->2. do the CSD search to extract cifs or just run existing files; for preexisting files set up the name for the topolgically labelled or raw cifs here
do_search = True

# ->3. set top_processed to True or False depending on whether you want to run john_top.py to do topology labelling
# ->3. 23/5/2020 changed things so top_processed means go through john_top.py
# ->3. john_top.py itself does various exclusions/group omissions.  relabelling the atoms is just one of the choices.  Set flags for this in #4
top_processed = True #go through john_top.py

# ->4. set up atoms/groups to replace by C then exclude here, applied in john_top and some (currently) to basicsearch if no topological analysis
# ->4. Sam rewrote ABn_a called sam_Abn_a in csdtools due to old script seeing any A connected to any 4B. Simplified FeCl5 for example and catena structures.
water_exclude = False # get rid of all water molecules - Sam now covered by small organics
replace_type_by_c = ['N', 'O'] # atoms that get changed to C because they'll be similar in size and could disorder
excluded_atoms = ['H','D'] # atoms to delete before symmetry analysis etc
ABn_list=['C_Cl_3', 'C_F_3', 'C_Br_3', 'Cl_O_4','Br_O_4','I_O_4', 'Cl_C_4','Br_C_4','I_C_4', 'B_F_4','B_H_4','Fe_Cl_4', 'Fe_Cl_5', 'Fe_Br_4','Re_O_4', 'Re_C_4', 'Ga_Cl_4', 'Mn_Br_4', 'Mn_Cl_4', 'Zn_Cl_4', 'Cu_Cl_4', 'Co_Cl_4', 'Al_Cl_4', 'Tl_Cl_4', 'Be_Cl_4', 'Cd_Cl_4', 'Zn_Br_4', 'Co_Br_4', 'Al_Br_4', 'Hg_Br_4', 'Ga_Br_4', 'Zn_I_4', 'Hg_I_4', 'In_I_4', 'Al_I_4'] # groups that likely to be disordered so B atoms removed - Sam -expanded by me, there could be more
topology_relabel = True # this chooses whether to use topology labelling or not

# ->5. set which entries in CSD to search and add to hit_list here; use 2000000 to search everything 0 to just use hit_list defined below.  Can search every nth entry if wanted
john_basicsearch.highest_entry = 0
john_basicsearch.nth_entry = 1

# ->6. pre-populate the hit_list here.  If the list is empty, just search database according to section 5 instructions.
# ->6. or can read from a .gcd file of refcodes; or use 32 known FE list; or a 311 test set or enter your own list
with open('known_ferros.gcd') as infile:
#with open('known_ferros.gcd') as infile:
    gcd_list = infile.readlines()
    gcd_list = [x.strip() for x in gcd_list]
john_basicsearch.hit_list = gcd_list
#john_basicsearch.hit_list = [] #comment out to include 32 knowns directly from john_basicsearch
#32 known FEs used for testing:
#john_basicsearch.hit_list = ['ABIQOU', 'ABIRAH', 'AMBACO07', 'ANPHPR01', 'BOLDIP11', 'CBUDCX', 'CYHDAC', 'FAWSAZ01', 'FIGNER', 'FINTAZ02', 'FINTIH', 'FOSGID', 'GLYCIN', 'GLYCIN01', 'GUMMUW', 'HORFAV', 'JAWQUT01', 'JAYPUU01', 'JIDBIK', 'JOJYOY','KOWYEA02', 'MAMPUM02', 'MAMPUM06', 'MEGJEP', 'REZBOP','REZBOP02', 'RIBLET', 'SIYWUT01', 'TEJKUO', 'TIVDAF01', 'UBOTIQ01', 'WIVDUA05']
john_basicsearch.hit_list = ['LIPZUI05']
print('Preset hit list is',john_basicsearch.hit_list)

# ->7. process pyroelectric pointgroups or all non-centrosymmetric: edit john_basicsearch.py around line 117

# ->8. make false to just do the search and skip running findsym
do_run_fs = True

# ->9. set findsym tolerances here; they default otherwise - Sam never messed with these
john_runfs.lat_tol      = 1.5    # latticeTolerance for cell parameters 1.5 optimal
john_runfs.atom_tol     = 1.8    # atomicPositionTolerance 1.8 optimal
john_runfs.atom_rel_tol = 0.999  # atomicPositionMaxTolerance atom tolerance relative to shortest A - A distance if atom_tol is too large i.e. min(atom_tol,atom_rel_tol*shortest_NN_dist) 0.99

# ->10. Set whether to interact with the database or not. Setting to false will still produce an excel file of results - Sam
interact_with_database = False
new_table = False   #if the interact is true and new table is false then the master table will be updated with new hits
#new_table_name = 'TestTable2'  # uncomment and name if new_table is true. make sure name doesn't already exist in the database or will crash
#john_runfs.new_table_name = new_table_name  # uncomment if new_table is true

# ->11. set whether to generate superimposed cifs of findsym output and input structures. - Sam
compare_cifs = True

if interact_with_database == True:
	run_comments = input('\nRun note --->')
	john_runfs.run_comments = run_comments


if interact_with_database == True and new_table == True:
	new_table_name = input('\nNew table name ---->')
	john_runfs.new_table_name = new_table_name
	con = sqlite3.connect('/home/dpqk92/PycharmProjects/csd_search_john/results_database.sqlite')
	cur = con.cursor()
	cur.execute("""CREATE TABLE {tab} (
				ID integer,
				Refcode text,
				LT_SpaceGroup text,
				HT_SpaceGroup text,
				AizuAllowed text,
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
# from here onwards run search, topology analysis (or just write raw cifs), runfs
# varius preferences from #1 to #9 above are sent to the other scripts
if delete_old_files:
	for item in glob.iglob('MasterOutputs/*/*'):
		os.remove(item)

john_basicsearch.excluded = excluded_atoms
john_basicsearch.replace_type_by_c = replace_type_by_c
john_top.excluded = excluded_atoms
john_top.replace_type_by_c = replace_type_by_c
john_top.ABn_list = ABn_list
john_top.water_exclude = water_exclude
john_top.topology_relabel = topology_relabel
john_runfs.new_table = new_table # Sam variables 
john_runfs.interact_with_database = interact_with_database
john_runfs.run_version = run_version

john_basicsearch.search()  # running the search and topological analysis

if top_processed:
	print("Running cifs through topological analysis")
	if do_search:
		john_basicsearch.run_top_an()
		john_basicsearch.write_CSD_cif() # writes out cifs as found in the database for use in john_compcif.py - sam
	john_runfs.cif_dir = sorted(glob.glob(top_out_dir+'/*.cif'))  # redefining some more directories
	john_runfs.top_processed = True
else:
	print("Will run raw cifs without topological analysis")
	if do_search:
		john_basicsearch.write_raw_cif() # option to write the raw cifs i.e. before topological analysis
	john_runfs.cif_dir = sorted(glob.glob(raw_out_dir+'/*.cif'))  # redefining some more directories
	john_runfs.top_processed = False

john_runfs.other_out_dir = other_out_dir
john_runfs.cif_out_dir = cif_out_dir
john_runfs.errant_out_dir = errant_out_dir
john_runfs.delete_unchanged_cif_inp_out = delete_unchanged_cif_inp_out
john_runfs.topology_relabel = topology_relabel
john_runfs.output_shifts = output_shifts

if do_run_fs:
	john_runfs.main()  # running the combined findsym script

if compare_cifs: #Sam - calls new py file for superimposing findsym structures onto raw cifs.
	with open('compcif_errors.txt', 'w+') as outfile:
		for file in os.listdir(cif_out_dir):
			try:
				filename = os.fsdecode(file)
				refcode = filename.strip('_top_findsym.cif')
				john_compcif3.main(refcode)
			except:
				outfile.write(str(refcode)+'\n')
				continue


end = time.time()

with open('MasterOutputs/Results.csv', 'a') as outfile:
	outfile.write('Script took '+str(end-start)+' seconds\n')

print('Script took '+str(end-start)+' seconds')


