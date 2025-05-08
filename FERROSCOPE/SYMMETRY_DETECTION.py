import sqlite3
import glob
import subprocess
import os
import re
import shutil
import ccdc.io
import time
import string
import numpy as np
import ast
import SUPERIMPOSE

# File set-up for running standalone. Otherwise, these are set by MASTER.py
# region
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

results = 'OutputCIFS/Results.csv'
results2 = 'OutputCIFS/Results2.csv'
errors = 'OutputCIFs/errors.txt'
notes_file = 'OutputCIFS/process_notes.txt'

lat_tol = 1.5  # lattice constant tolerance
atom_tol = 1.8  # atom position tolerance
atom_rel_tol = 0.999  # atom tolerance relative to shortest A - A distance if atom_tol is too large i.e. min(atom_tol,atom_rel_tol*shortest_NN_dist)
delete_unchanged_cif_inp_out = False

# endregion

# Symmetry dictionaries
# region
pg_spg = {
	'pg_1': ['1'],
	'pg_1bar': ['2'],
	'pg_2': ['3', '4', '5'],
	'pg_m': ['6', '7', '8', '9'],
	'pg_2upm': ['10', '11', '12', '13', '14', '15'],
	'pg_222': ['16', '17', '18', '19', '20', '21', '22', '23', '24'],
	'pg_mm2': ['25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46'],
	'pg_mmm' : ['47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74'],
	'pg_4' : ['75', '76', '77', '78', '79', '80'],
	'pg_4bar' : ['81', '82'],
	'pg_4upm' : ['83', '84', '85', '86', '87', '88'],
	'pg_422' : ['89', '90', '91', '92', '93', '94', '95', '96', '97', '98'],
	'pg_4mm' : ['99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110'],
	'pg_4bar2m' : ['111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122'],
	'pg_4upmmm' : ['123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142'],
	'pg_3' : ['143', '144', '145', '146'],
	'pg_3bar' : ['147', '148'],
	'pg_32' : ['149', '150', '151', '152', '153', '154', '155'],
	'pg_3m' : ['156', '157', '158', '159', '160', '161'],
	'pg_3barm' : ['162', '163', '164', '165', '166', '167'],
	'pg_6' : ['168', '169', '170', '171', '172', '173'],
	'pg_6bar' : ['174'],
	'pg_6upm' : ['175', '176'],
	'pg_622' : ['177', '178', '179', '180', '181', '182'],
	'pg_6mm' : ['183', '184', '185', '186'],
	'pg_6barm2' : ['187', '188', '189', '190'],
	'pg_6upmmm' : ['191', '192', '193', '194'],
	'pg_23' : ['195', '196', '197', '198', '199'],
	'pg_m3bar' : ['200', '201', '202', '203', '204', '205', '206'],
	'pg_432' : ['207', '208', '209', '210', '211', '212', '213', '214'],
	'pg_4bar3m' : ['215', '216', '217', '218', '219', '220'],
	'pg_m3barm' : ['221', '222', '223', '224', '225', '226', '227', '228', '229', '230']
}

Aizu_transitions = {
	'pg_1': ['pg_1bar', 'pg_2', 'pg_m', 'pg_2upm', 'pg_222', 'pg_mm2', 'pg_mmm', 'pg_4', 'pg_4bar', 'pg_4upm', 'pg_422',  'pg_4mm', 'pg_4bar2m', 'pg_4upmmm', 'pg_3', 'pg_3bar', 'pg_32', 'pg_3m', 'pg_3barm', 'pg_6', 'pg_6bar', 'pg_6upm', 'pg_622', 'pg_6mm', 'pg_6barm2', 'pg_6upmmm', 'pg_23', 'pg_m3bar', 'pg_432', 'pg_4bar3m', 'pg_m3barm'],
	'pg_2': ['pg_2upm', 'pg_222', 'pg_4bar', 'pg_422', 'pg_4bar2m', 'pg_32', 'pg_3barm', 'pg_622', 'pg_23', 'pg_432'],
	'pg_m': ['pg_2upm', 'pg_mm2', 'pg_mmm', 'pg_4upm', 'pg_4mm', 'pg_4bar2m', 'pg_4upmmm', 'pg_3m', 'pg_3barm', 'pg_6bar', 'pg_6upm', 'pg_6mm', 'pg_6barm2', 'pg_6upmmm', 'pg_4bar3m', 'pg_m3barm'],
	'pg_mm2': ['pg_mmm', 'pg_4bar2m', 'pg_4upmmm', 'pg_6barm2', 'pg_6upmmm', 'pg_m3bar', 'pg_4bar3m', 'pg_m3barm'],
	'pg_4': ['pg_4upm', 'pg_422', 'pg_432'],
	'pg_4mm': ['pg_4upmmm', 'pg_m3barm'],
	'pg_3': ['pg_3bar', 'pg_32', 'pg_6bar', 'pg_23', 'pg_432'],
	'pg_6': ['pg_6upm', 'pg_622'],
	'pg_3m': ['pg_3barm', 'pg_6barm2', 'pg_4bar3m', 'pg_m3barm'],
	'pg_6mm': ['pg_6upmmm']
}

space_group_centering = {
	'P' : ['1', '2', '3', '4', '6', '7', '10', '11', '14', '16', '17', '18', '19', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '75', '76', '77', '78', '81', '82', '83', '84', '85', '86', '89', '90', '91', '92', '93', '94', '95', '96', '99', '100', '101', '102', '103', '104', '105', '106', '111', '112', '113', '114', '115', '116', '117', '118', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '143', '144', '145', '147', '149', '150', '151', '152', '153', '154', '156', '157', '158', '159', '162', '163', '164', '165', '168', '169', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', '191', '192', '193', '194', '195', '198', '200', '201', '205', '207', '208', '212', '213', '215', '218', '221', '222', '223', '224'],
	'C' : ['5', '8', '9', '12', '13', '15', '20', '21', '35', '36', '37', '38', '39', '40', '41', '63', '64', '65', '66', '67', '68'],
	'F' : ['22', '42', '43', '69', '70', '196', '202', '203', '209', '210', '216', '219', '225', '226', '227', '228'],
	'I' : ['23', '24', '44', '45', '46', '71', '72', '73', '74', '79', '80', '87', '88', '97', '98', '107', '108', '109', '110', '119', '120', '121', '122', '139', '140', '141', '142', '197', '199', '204', '206', '211', '214', '217', '220', '229', '230'],
	'R' : ['146', '148', '155', '160', '161', '166', '167'],
}

primitive_elements = {
	'1' : ['pg_1'],
	'2' : ['pg_1bar', 'pg_2', 'pg_m'],
	'3' : ['pg_3'],
	'4' : ['pg_2upm', 'pg_222', 'pg_mm2', 'pg_4', 'pg_4bar'],
	'6' : ['pg_3bar', 'pg_32', 'pg_3m', 'pg_6', 'pg_6bar'],
	'8' : ['pg_mmm', 'pg_4upm', 'pg_422', 'pg_4mm', 'pg_4bar2m'],
	'12': ['pg_3barm', 'pg_6upm', 'pg_622', 'pg_6mm', 'pg_6barm2', 'pg_23'],
	'16': ['pg_4upmmm'],
	'24': ['pg_6upmmm', 'pg_m3bar', 'pg_432', 'pg_4bar3m'],
	'48': ['pg_m3barm']
}


#point group categories
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
print('Number of centrosymmetric space groups len(centro_spacegroups)= ',len(centro_spacegroups))
print('Number of non centrosymmetric space groups len(non_centro_spacegroups) =',len(non_centro_spacegroups))
print('Number of non pyrolelectric space groups len(pyropgs) =',len(pyropgs))

# endregion


class IndexTooHigh(Exception):
	"""A custom exception, raised when our symmetry index is too large"""
	pass


class NoAtomsError(Exception):
	"""A custom exception, raised when our we have excluded every atom in the CIF"""
	pass


# an Atom class for get_shifts
class Atom:
	def __init__(self):
		self.name = ''
		self.type = ''
		self.oldtype = ''
		self.wyckoff =''
		self.id = 0
		self.x = 0.0
		self.y = 0.0
		self.z = 0.0
		self.Z = 0
		self.mult = 0
		self.occ = 1
		self.rms = 0
		self.max = 0


# get the rms shifts for atoms from findsym output and corresponding old label types
def get_shifts(findsym_out,shift,old_labels,out_shift):
	max_rms=0
	rmstotal=0
	unique=0
	global atoms
	atoms = []
	maxtype=''
	maxoldtype=''
	for i,line in enumerate(findsym_out):
		if re.search('Wyckoff ',line):
			firstline=line
			coords=findsym_out[i + 1].split()
			startline=i
			unique +=1
		if re.search('RMS position',line):
			atom=Atom()
			atom.mult=float(i - (startline + 1))
			atom.x=float(coords[1])
			atom.y=float(coords[2])
			atom.z=float(coords[3])
			atom.name=re.search(r'\((.*?)\)',firstline).group(1)
			atom.type=atom.name[0:3].upper()
			try:
				atom.oldtype=old_labels[atom.type]
			except:
				atom.oldtype='UN'
			atom.wykcoff=firstline.split()[2]
			atom.rms=float((line.split()[4]).replace(',',''))
			atom.max=float(line.split()[9])
			atoms.append(atom)
		if re.search('Ratio ', line):
			R = int(findsym_out[i].split()[-1])

	sum_mi = 0
	sum_mi_rms2 = 0
	for atom in atoms:
		if atom.max > max_rms:
			max_rms = atom.max
			maxtype = atom.type
			maxoldtype = atom.oldtype
		sum_mi = sum_mi+atom.mult
		sum_mi_rms2 = sum_mi_rms2 + atom.mult * atom.rms**2
		rmstotal=rmstotal+atom.rms
	average_rms = rmstotal/len(atoms)
	overall_rms = (sum_mi_rms2/sum_mi)**0.5
	overall_rss_C = (R*sum_mi_rms2)**0.5
	overall_rss_P = overall_rss_C / R ** 0.5

	if out_shift:
		with open(shift,'w') as outfile:
			outfile.write('name\ttype\toriginal_type\tmult\tx\ty\tz\trms_shift\tmax_rms_shift\n')
			for atom in atoms:
				outfile.write("%-4s\t%-4s\t%-4s\t%3.1f\t%10f\t%10f\t%10f\t%10f\t%10f\n" % (atom.name,atom.type,atom.oldtype,atom.mult,atom.x,atom.y,atom.z,atom.rms,atom.max))
	return overall_rms,max_rms,maxtype,maxoldtype, overall_rss_C, overall_rss_P


def space_group_finder(spg):
	#Sam function for finding point group from space group using new dictionaries at top of file
	for x,y in pg_spg.items():
		if spg in y:
			return(x)
		else:
			pass
	else:
		print("Couldn't find point group from space group")


def aizu_allowed(old, new):
	#Sam function to check if the transition found by FINDSYM is a full Aizu para-ferro transition
	potentials = Aizu_transitions.get(old)
	if potentials != None: #if loop included to prevent crash due to potentials being empty after space_group_finder fails
		if new in potentials:
			return True
		else:
			return False
	else:
		print('Not adding to list')
		return False


def clean_up(cif):
	"""This function replaces the atomic types with labels given by topological analysis and then sorts the CIF by label
	Returns a dictionary with the new labels as keys and the old as entries, also the key 'last_index' gives the number of
	symmetry indices in the CIF"""
	alphabet = string.ascii_uppercase
	with open(cif, 'r') as infile:  # opens the cif and puts all the lines in a list
		total = infile.readlines()
	try:
		total.remove('_atom_site_U_iso_or_equiv\n')     #sam added these two lines after 2022 upgrade. these weren't in cif files before upgrade I think
		total.remove('_atom_site_thermal_displace_type\n')
	except ValueError:
		print("No thermal displacement lines")
		pass
	cut2 = []  # must define it first in case something goes wrong
	for i, line in enumerate(total):
		if line == '_atom_site_fract_z\n':  # finds the first line of atomic coordinates
			cut1 = total[i + 1:]
			found = False
			start = i + 1
			for j, line2 in enumerate(cut1):  # finds the last line and cuts just the required section
				if line2 == 'loop_\n' and found == False:
					cut2 = cut1[:j]
					found = True
				elif line2 == '\n' and found == False:
					cut2 = cut1[:j]
					found = True

	if cut2 == []:  # Just in case we somehow remove all the atoms in the structure
		print('No atoms present')
		raise NoAtomsError

	for item in cut2:  # here we remove the old coordinates from the cif before we add the new ones, must be before we change atom type
		total.remove(item)  # total is the list of the whole cif

	new_coords = np.chararray((len(cut2), 5), itemsize=50, unicode=True)  # creates a char array so we can sort the coordinates need unicode=True for python3.7
	corresponding_label = dict([])  # creates a dictionary for remembering the old atom types
	for i, line in enumerate(cut2):  # this section finds the sym index and converts it to the right letter code
		separated = line.split()
		separated = separated[0:5] #sam included this to remove temperature factors from cif file
		label = separated[0]
		index = int(label.split('_')[1])
		old_type = separated[1]
		v = index - 1  # since our indexes start at 1 but python starts from 0
		max_power = 3  # No. of characters in the new label, the max number of groups is (26**max_power), 3 is 17576
		power = max_power
		new_type = ''

		if index > 26 ** max_power:  # accounts for having too many symmetry groups
			raise IndexTooHigh

		while v > 0:  # finds the correct letter code for each number
			if v < 26**power:
				power -= 1
			else:
				this_letter = v//(26**power)  # // is integer division eg. 51//26 = 1
				new_type += alphabet[this_letter]
				v -= this_letter * (26**power)

		new_type += 'A' * power  # accounts for multiples of 26 directly
		new_type = 'A' * (max_power - len(new_type)) + new_type  # Adds any leading A's to make the right length code

		separated[1] = new_type  # changes the atom type to the new letter
		corresponding_label[new_type] = old_type
		for j in range(len(separated)):
			new_coords[i, j] = separated[j]  # adds the whole line to the array

	new_coords = new_coords[new_coords[:, 1].argsort()]  # sorts the coords by new type
	last_index = new_coords[-1].split()[0][0].split('_')[1]  # finds the highest number index for the results file
	corresponding_label['last_index'] = last_index

	sorted = []
	for i in range(len(cut2)):  # reconverts the array to a standard list
		temp = ' '.join(new_coords[i, :])
		sorted.append(temp + '\n')

	for j, item in enumerate(sorted):  # then we add the new coords to the cif file
		total.insert(start + j, item)

	with open(cif, 'w') as outfile:  # and finally write out the changes
		outfile.writelines(total)

	return corresponding_label  # returns the dictionary of new/old atom types for delabelling later


def delabel(cif, dictionary, cif_type=2):
	"""A function to convert the atom types in a CIF back to the original elements from the symmetry indices
	dictionary contains the new atom type as a key for the old atom type as entries.
	Type=1 delabels cifs from CSD, type 2 delabels cifs from findsym"""

	with open(cif, 'r') as infile:  # opens the cif and puts all the lines in a list
		total = infile.readlines()

	cut2 = []
	for i, line in enumerate(total):
		if cif_type == 1:
			if line == '_atom_site_fract_z\n':  # finds the first line of atomic coordinates
				cut1 = total[i + 1:]
				found = False
				start = i + 1
				for j, line in enumerate(cut1):  # finds the last line and cuts just the required section
					if line == 'loop_\n' and found == False:
						cut2 = cut1[:j]
						found = True
						end = start + j
					elif line == '\n' and found == False:
						cut2 = cut1[:j]
						found = True

		elif cif_type == 2:
			if line == '_atom_site_fract_symmform\n':  # finds the first line of atomic coordinates
				cut1 = total[i + 1:]
				found = False
				start = i + 1
				for j, line in enumerate(cut1):  # finds the last line and cuts just the required section
					if line == 'loop_\n' and found == False:
						cut2 = cut1[:j]
						found = True
						end = start + j
					elif line == '\n' and found == False:
						cut2 = cut1[:j]
						found = True

	new_str = []
	for i, line in enumerate(cut2):
		separated = line.split()
		new_type = separated[1].upper()  # finds the new atomic type and switches it for the old one, need upper() as fs converts to lower
		separated[1] = dictionary[new_type]
		new_str.append(' '.join(separated) + '\n')

	for line in cut2:  # here we remove all the old coords from the cif file
		total.remove(line)

	for j, item in enumerate(new_str):  # then we add the new coords to the cif file
		total.insert(start + j, item)

	with open(cif, 'w') as outfile:  # and finally write out the changes
		outfile.writelines(total)


def main():
	"""This function firstly runs all the cif files through findsym_cifinput and then through findsym itself
		27/4/2020 JSOE added top_processed depending on whether raw cifs or not """
	with open(results, 'w+') as outfile:
		outfile.write('ref_code\told_sym_no\tnew_sym_no\tname\tformula\tno_of_atoms\tno_of_symmetry_groups\ta_in\tb_in\tc_in\tal_in\tbe_in\tga_in\t'
					  'vol_in\ta_out\tb_out\tc_out\tal_out\tbe_out\tga_out\tvol_out\tauthor(s)\tjournal_name\tjournal_volume\tpub_year'
					  '\tpage_no\tDOI\n')
		outfile.close()
	with open(results2, 'w+') as outfile:
		outfile.write('serial\tref_code\thyperlink\told_sym_no\tnew_sym_no\told_sgname\tnew_sgname\tAizu_allowed\tcentro\toverall_rms\tmax_rms\tmax_type\tmax_oldtype\tvol_ratio\tname\tformula\tno_of_atoms\tno_of_symmetry_groups\ta_in\tb_in\tc_in\tal_in\tbe_in\tga_in\t'
					  'vol_in\ta_out\tb_out\tc_out\tal_out\tbe_out\tga_out\tvol_out\tauthor(s)\tjournal_name\tjournal_volume\tpub_year'
					  '\tpage_no\tDOI\tT_datacoll\tT_melt\tprocess_note\n')
		outfile.close()
	with open(errors, 'w+') as outfile:
		outfile.write('Sample_name\tError\n')
		outfile.close()

	attempted = 0  # counter to check how many samples it has run through
	failed = 0  # simple counter to check how many input files fail in some way
	changed = 0  # counter to check how many have a higher symmetry form found
	changed_centro = 0 # counter for number of changed which are to centrosymmetric
	hitlist = []  # list of refcodes that are changed
	crashlist = []  # list of refcodes that have crashed
	entries_updated = 0 #counter for entries in database updated
	entries_added = 0 #counter for entries in database added
# read in a dictionary with process notes for each entry, send a different default message if file missing
	notes_dict = {}
	try:
		with open(notes_file,'r') as f:
			default='No note found for entry'
			for line in f:
				(key, val) = line.split(':')
				notes_dict[key] = val.strip()
	except OSError:
		default='No file with notes read'
	for file in cif_dir:
		crash = False
		end_line = False
		attempted += 1
		cif = os.path.realpath(file)  # finds the absolute path to the input cif
		file_name = os.path.basename(file)  # takes just the filename from the whole path
		if top_processed: #this is still top_processed not top_relabel for now 22/5/2020
			ref_code = file_name[0:-8]  # defines the refcode of the CSD entry i.e. it removes the '_top.cif'
		else:
			ref_code = file_name[0:-4]  # defines the refcode of the CSD entry i.e. it removes the '_top.cif'
		inp = file_name.replace('.cif', '.inp')
		inp = os.path.join(other_out_dir, inp)  # creates a new string with the same file name but the .inp ext
		out = file_name.replace('.cif', '.out')
		out = os.path.join(other_out_dir, out)  # creates a new string with the same file name but the .out ext
		shift = file_name.replace('.cif', '.shift')
		shift = os.path.join(other_out_dir,shift)
		findsym_cif = file_name.replace('.cif', '_findsym.cif')
		findsym_cif = os.path.join(cif_out_dir, findsym_cif)  # creates a new file name for the findsym_cif output

		if topology_relabel:
			print('\nRelabelling CIF file {} of {} to AAA type names ready for findsym'.format(attempted,len(cif_dir)))
			try:
				old_labels = clean_up(file)  # relabels and sorts the atoms
				print('Old labels dictionary is',old_labels)

			except IndexTooHigh:  # exception for when we have more than 702 symmetry groups
				failed += 1
				crashlist.append(ref_code)
				err_line = 'Error: Too many symmetry groups in unit cell'
				print(err_line+'\n')
				with open(errors, 'a') as outfile:  # writes info to the error file
					outfile.write(ref_code+'\t'+err_line+'\n')
					outfile.close()
				shutil.copy(cif, errant_out_dir)  # copies the top_an cif to a new directory
				continue  # skips the rest of the code bc it won't function

			except NoAtomsError:  # exception for when we have excluded all of the atoms in the structure
				failed += 1
				crashlist.append(ref_code)
				err_line = 'Error: No atoms left in the structure after exclusions'
				print(err_line+'\n')
				with open(errors, 'a') as outfile:  # writes info to the error file
					outfile.write(ref_code+'\t'+err_line+'\n')
					outfile.close()
				shutil.copy(cif, errant_out_dir)  # copies the top_an cif to a new directory
				continue  # skips the rest of the code bc it won't function
		else:
			old_labels = []

		with open(inp, 'w+') as outfile:  # this block runs the cif into 'findsym_cifinput' to get the .inp file
			subprocess.call(['findsym_cifinput', cif], stdout=outfile)
			outfile.close()
			print('Findsym input file created from ' + cif)

		# reports input spacegroup, try a few different labels in the cif to maximise chance of success
		with open(cif, 'r') as infile:
			for line in infile:
				if re.search('sym.*number', line):
					print('Input symmetry:' + line.rstrip())
					old_sym_no = int(re.findall(r'\d+', line)[0])  # finds the space group number from the input cif
				if re.search('space.*number', line):
					print('Input symmetry:' + line.rstrip())
					old_sym_no = int(re.findall(r'\d+', line)[0])

		with open(inp, 'r') as infile:
			text = infile.read()
			# text = re.sub('^#.*\n','', text)
			text = re.sub(re.compile('^#.*\n', re.MULTILINE), '', text)
			text = re.sub('!latticeTolerance.*\n.*', '!latticeTolerance\n %f' % lat_tol, text)
			text = re.sub('!atomicPositionTolerance.*\n.*', '!atomicPositionTolerance\n %f' % atom_tol, text)
			text = re.sub('!atomicPositionMaxTolerance.*\n.*', '!atomicPositionMaxTolerance\n %f' % atom_rel_tol, text)
			infile.close()
		with open(inp, 'w') as outfile:
			outfile.write(text)

		with open(out, 'w+') as outfile:  # running findsym with the new timer
			task = subprocess.Popen(['findsym', inp], stdout=outfile)  # starts the process
			timer = 60.0
			while timer > 0 and task.poll() is None:  # a primitive timer, stops when t hits 0 or we get a result
				time.sleep(0.5)
				timer -= 0.5

			if timer > 0 and task.poll() is not None:  # here we get a result before timer runs out
				crash = False
			elif timer <= 0 and task.poll() is None:  # here the timer runs out before a result - taken as an error
				crash = True
				err_line = 'This file timed out'
			else:  # in theory this should never happen but it's included just in case
				crash = True
				err_line = 'I don\'t know why this happened - error with the timer'

		with open(out, 'r') as outfile:  # report output spacegroup and check for crashes
			for line in outfile:
				if re.search('Space Group', line):
					print('Output symmetry:                   ' + line.rstrip())
					new_sym_no = int(re.findall(r'\d+', line)[0])  # finds the space group from the new cif
					crash = False
					break
				if re.search('Error ', line):
					crash = True
					err_line = line.rstrip()
					break
				if re.search('ombed', line):
					crash = True
					err_line = line.rstrip()
					break
				# else:  # if we do not hit any of the above it's likely the file timed out
				#     crash = True
				#     err_line = 'An unknown error occurred'
			if crash == True:
				print('**Check .inp and .out files carefully: findsym may have crashed**')
				failed += 1
			if crash == False:
				print('File put through findsym successfully')

# put johns extra loop here to check for no "end of cif" if its not already flagged as a crash
		if crash == False:
			with open(out,'r') as outfile:
				if not 'end of cif' in outfile.read():
					#print('Found end of in file')
				#else:
					print('....but didnt find "end of cif" in the output file')
					crash = True
					failed +=1
					err_line = 'No error code but didnt find "end of cif" in the output file'

#
		if crash == True:  # writes the error and then moves to the next entry
			crashlist.append(ref_code)
			print(err_line+'\n')
			with open(errors, 'a') as outfile:
				outfile.write(ref_code+'\t'+err_line+'\n')
				outfile.close()
			shutil.copy(cif, errant_out_dir)  # copies the top_an cif and inp/out files to a new directory
			shutil.copy(inp, errant_out_dir)
			shutil.copy(out, errant_out_dir)
			continue

		# make a clean output cif file
		tempcif = []
		with open(out, 'r') as infile:
			findsym_out = infile.readlines() #john changed 27/5 use to just use infile
		copy = False
		for line in findsym_out:
			if re.search('CIF file', line):
				copy = True
			if copy == True:
				tempcif.append(line.rstrip())
		#section here to get rms shifts etc from out
		#call shifts here
		#get_shifts(findsym_out,ref_code,old_labels,output_shifts)
		overall_rms, max_rms, max_type, max_oldtype, overall_rss_C, overall_rss_P = get_shifts(findsym_out,shift,old_labels,output_shifts)

		with open('rss.txt', 'a+') as outfile:
			outfile.write(ref_code + ' '+ str(overall_rss_C) + ' ' + str(overall_rss_P) + '\n')

		#section here to check for tolerances hit
		if topology_relabel:
			dashed_lines = []
			dashed_lines_num = 0
			line_num = 0
			tolerance_output = ""
			for line in findsym_out:
				if re.search('------------------------------------------', line):
					dashed_lines_num += 1
					dashed_lines.append(line_num)
				line_num+=1
			if dashed_lines_num == 3: #we've hit tolerances
				for num in range(dashed_lines[0], dashed_lines[1]):
					if findsym_out[num][0] == "A":
						tolerance_line = findsym_out[num].strip()
						new_label = old_labels[tolerance_line[0:3]]
						tolerance_line = tolerance_line + " " + new_label
						tolerance_output = tolerance_output + tolerance_line + '\n'
			elif dashed_lines_num == 2: #we've not hit tolerances
				pass
			else:
				with open('errors.txt', 'a+') as outfile:
					outfile.write(str(molecule.identifier) + ' : SYMMETRY_DETECTION.tolerance_hitting: ' + str(e) + '\n')
			print(tolerance_output)
		else:
			dashed_lines = []
			dashed_lines_num = 0
			line_num = 0
			tolerance_output = ""
			for line in findsym_out:
				if re.search('------------------------------------------', line):
					dashed_lines_num += 1
					dashed_lines.append(line_num)
				line_num += 1
			if dashed_lines_num == 3:  # we've hit tolerances
				for num in range(dashed_lines[0], dashed_lines[1]):
					if findsym_out[num][0] == "t":
						type_line = num
				for num2 in range(type_line + 1, dashed_lines[1]):
					tolerance_line = findsym_out[num2].strip()
					new_label = tolerance_line[0:3]
					tolerance_line = tolerance_line + " " + new_label
					tolerance_output = tolerance_output + tolerance_line + '\n'
			elif dashed_lines_num == 2:  # we've not hit tolerances
				pass
			else:
				with open('errors.txt', 'a+') as outfile:
					outfile.write(
						str(molecule.identifier) + ' : SYMMETRY_DETECTION.tolerance_hitting: ' + str(e) + '\n')
			print(tolerance_output)

		with open(out, 'r') as infile:  # report output spacegroup and check for crashes
			findsym_out = infile.readlines()
		with open(findsym_cif, 'w+') as outfile:
			for item in tempcif:
				outfile.write('%s\n' % item)
		if topology_relabel:
			try:
				delabel(cif, old_labels, 1)  # reverts the original and findsym cifs back to the proper atom types
				delabel(findsym_cif, old_labels, 2)

			except KeyError:  # when there is a problem with the key it is likely because of an error in findsym
				failed += 1
				crashlist.append(ref_code)
				err_line = 'Error in delabelling the CIF\n'
				print(err_line)
				with open(errors, 'a') as outfile:
					outfile.write(ref_code + '\t' + err_line)
				shutil.copy(cif, errant_out_dir)  # copies the top_an cif and inp/out files to a new directory
				shutil.copy(inp, errant_out_dir)
				shutil.copy(out, errant_out_dir)
				continue

		with open(cif, 'r') as infile:  # finds the old cell params
			text_temp = infile.readlines()
			for i, line in enumerate(text_temp):
				if re.search('_cell_volume', line):
					vol_in = line.split()[1]
					vol_in = re.sub(r'\(.*?\)', "", vol_in)

		with open(findsym_cif, 'r') as infile:  # finds the new cell params
			text_temp2 = infile.readlines()
			for i, line in enumerate(text_temp2):
				if re.search('_cell_volume', line):
					vol_out = line.split()[1]

		#check if volume changes so spg --> spg transitions don't get lost
		R_centred = [146, 148, 155, 160, 161, 166, 167]

		vol_ratio = float(vol_in)/float(vol_out)
		vol_100 = vol_ratio * 100
		vol_round = int(round(vol_100,))
		if old_sym_no in R_centred:
			if old_sym_no == new_sym_no and str(vol_ratio)[:4] == "0.33":
				vol_round = 100  #findsym changing cubic axes to hex axes counted as a volume change and thus a hit, this fixes



		if old_sym_no != new_sym_no or vol_round not in range(95, 105):  # adds entries that change to a results file
			print('This crystal has a higher symmetry form found by findsym')
			print('Overall_rms_shift:',overall_rms,'Max_rms_shift:',max_rms,'Type:',max_type,'Original_type:',max_oldtype)
			aizu_yes = ''
			old_spg = str(old_sym_no)
			new_spg = str(new_sym_no)
			old_pg = space_group_finder(old_spg)
			new_pg = space_group_finder(new_spg)
			if aizu_allowed(old_pg, new_pg) == True: #Sam - checking aizu allowed
				print('Aizu allowed - adding to results file')
				aizu_yes = 'Yes'
			else:
				print('Aizu disallowed - not adding to results file')
				aizu_yes = 'No'
				#continue
			if str(new_sym_no) in centro_spacegroups:
				centro = 'Centro'
				print('Higher symmetry form is centrosymmetric')
				changed_centro +=1
			else:
				print('Higher symmetry form is noncentrosymmetric')
				centro = 'Noncentro'
			changed += 1
			CSD = ccdc.io.EntryReader('csd')  # now we collect some information on the entry to write to the results
			hitlist.append(ref_code)
			entry = CSD.entry(ref_code)
			crystal = CSD.crystal(ref_code)
			name = entry.chemical_name
			temperature =entry.temperature
			if not temperature:
				temperature = 'Unknown Tcol'
			melting = entry.melting_point
			if not melting:
				melting ='Unknown Tmelt'
			#print('Refcode',ref_code,'at temperature',temperature)
			#print('Refcode',ref_code,'melting point',melting)
			if name != None:
				#john commented this out as didnt see why it was byte encoded, added pass
				#name = name.encode('utf8')
				pass
			formula = entry.formula
			if formula != None:
				#formula = formula.encode('utf8')
				pass
			author = entry.publication[0]
			author = str(author)
			if author != None:
				#author = author.encode('utf8')
				pass
			# collect data from CSD entry
			# region
			journal_name = entry.publication[1] #Sam added more of these to fill out the database
			journal_vol = entry.publication[2]
			year = entry.publication[3]
			page_no = entry.publication[4]
			doi = str(entry.publication[5])
			Z = crystal.z_value
			Zprime = crystal.z_prime
			calc_density = crystal.calculated_density
			analogue = entry.analogue
			bioactivity = entry.bioactivity
			ccdc_number = entry.ccdc_number
			colour = entry.color
			database = entry.database_name
			depositiondate = entry.deposition_date
			disorder = entry.disorder_details
			habit = entry.habit
			organic = entry.is_organic
			polymeric = entry.is_polymeric
			powder = entry.is_powder_study
			phase_transition = entry.phase_transition
			polymorph = entry.polymorph
			pressure = entry.pressure
			r_factor = entry.r_factor
			radiation = entry.radiation_source
			remarks = entry.remarks
			solvent = entry.solvent
			source = entry.source
			synonyms = entry.synonyms
			#endregion
			if topology_relabel:
				sym_group_no = old_labels['last_index']
			else:
				sym_group_no = 0
			with open(cif, 'r') as infile:  # finds the old cell params
				text_temp = infile.readlines()
				for i, line in enumerate(text_temp):
					if re.search('_cell_length_a', line):
						a_in = line.split()[1]
						a_in = re.sub(r'\(.*?\)',"",a_in)
					elif re.search('_cell_length_b', line):
						b_in = line.split()[1]
						b_in = re.sub(r'\(.*?\)',"",b_in)
					elif re.search('_cell_length_c', line):
						c_in = line.split()[1]
						c_in = re.sub(r'\(.*?\)',"",c_in)
					elif re.search('_cell_angle_alpha', line):
						al_in = line.split()[1]
						al_in = re.sub(r'\(.*?\)',"",al_in)
					elif re.search('_cell_angle_beta', line):
						be_in = line.split()[1]
						be_in = re.sub(r'\(.*?\)',"",be_in)
					elif re.search('_cell_angle_gamma', line):
						ga_in = line.split()[1]
						ga_in = re.sub(r'\(.*?\)',"",ga_in)
					elif re.search('_cell_volume', line):
						vol_in = line.split()[1]
						vol_in = re.sub(r'\(.*?\)',"",vol_in)
					elif re.search('_symmetry_space_group_name_H-M', line):
						old_sgname='???'
						#need the if as there is one entry in database (SACNIC) without proper space group
						if re.findall(r"'(.*?)'",line):
							old_sgname = re.findall(r"'(.*?)'",line)[0]
							old_sgname = re.sub(r"\s","_",old_sgname)
						print('Original HM spacegroup is:',old_sgname)


			with open(findsym_cif, 'r') as infile:  # finds the new cell params
				text_temp2 = infile.readlines()
				for i, line in enumerate(text_temp2):
					if re.search('_cell_length_a', line):
						a_out = line.split()[1]
					elif re.search('_cell_length_b', line):
						b_out = line.split()[1]
					elif re.search('_cell_length_c', line):
						c_out = line.split()[1]
					elif re.search('_cell_angle_alpha', line):
						al_out = line.split()[1]
					elif re.search('_cell_angle_beta', line):
						be_out = line.split()[1]
					elif re.search('_cell_angle_gamma', line):
						ga_out = line.split()[1]
					elif re.search('_cell_volume', line):
						vol_out = line.split()[1]
					elif re.search('_symmetry_space_group_name_H-M', line):
						new_sgname = re.findall(r'"(.*?)"',line)[0]
						new_sgname = re.sub(r"\s","_",new_sgname)
						print('Final HM spacegroup is:   ',new_sgname)

			with open(inp, 'r') as infile:  # finds the number of atoms in the unit cell
				text_temp3 = infile.readlines()
				for i, line in enumerate(text_temp3):
					if re.search('!atomCount', line):
						atom_no = text_temp3[i+1].strip()
			note=notes_dict.get(ref_code, default)
			print('Processing note is:',note)
			ref_code_link = 'https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid={0}&DatabaseToSearch=Published'.format(ref_code)
			doi_link = '=HYPERLINK("https://dx.doi.org/{}","DOI")'.format(doi)  # these format the refcode/DOI as clickable hyperlinks in excel

			with open(results, 'a') as outfile:  # write all the information to the results file
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}'
							  '\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n'.format(ref_code_link,
								old_sym_no, new_sym_no, name, formula, atom_no, sym_group_no, a_in, b_in, c_in, al_in, be_in, ga_in, vol_in,
								a_out, b_out, c_out, al_out, be_out, ga_out, vol_out, author, journal_name, journal_vol, year, page_no, doi_link))
			with open(results2, 'a') as outfile:  # write all the information to the results file
				outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:3.6f}\t{:3.6f}\t'
								'{}\t{}\t{:3.1f}\t{}\t{}\t{}\t{}\t'
								'{}\t{}\t{}\t{}\t{}\t{}\t{}\t'
								'{}\t{}\t{}\t{}\t{}\t{}\t{}\t'
								'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
								.format(changed, ref_code, ref_code_link, old_sym_no, new_sym_no, old_sgname, new_sgname, aizu_yes, centro, overall_rms, max_rms,
								max_type, max_oldtype, float(vol_in)/float(vol_out), name, formula, atom_no, sym_group_no,
								a_in, b_in, c_in, al_in, be_in, ga_in, vol_in,
								a_out, b_out, c_out, al_out, be_out, ga_out, vol_out,
								author, journal_name, journal_vol, year, page_no, doi_link,temperature,melting,note))

			# Sam's script for checking whether continuous and its respective irrep

			#region
			irrep = '-'
			continuous = False
			continuous_conditions = ['RG', 'yes']
			iso_input = 'VALUE PARENT {}\nVALUE SUBGROUP {}\nSHOW IRREP\nSHOW CONTINUOUS\nSHOW BASIS\nDISPLAY ISOTROPY\n'.format(new_sym_no, old_sym_no)
			iso_input = iso_input.encode()

			with open('isotropyout.txt', 'w+') as outfile:
				subprocess.run(['iso'], input=iso_input, stdout=outfile, stderr=subprocess.DEVNULL)
				outfile.close()

			with open('isotropyout.txt', 'r+') as infile:
				lines = infile.readlines()
				important_lines = lines[5:(len(lines) - 1)]
				infile.close()
			for line in important_lines:
				irrep = line.split()[0]
				continuous_allowed = line.split()[1]
				matrix = (line.split()[2])
				try:
					array = np.array(ast.literal_eval(matrix))
					zeros = (np.count_nonzero(array == 0))
					ones = (np.count_nonzero(array == 1) + np.count_nonzero(array == -1))
					if zeros == 6 and ones == 3 and continuous_allowed in continuous_conditions:  # 2nd order FE phast transition needs to be at the GM point and so the matrix must only three 1's
						continuous = True
						break
					else:
						irrep = '-'
						continuous = False
				except ValueError as error: # matrix is not purely 1's and 0's and is thus not continuous (in a ferroelectric sense anyway)
					continue
			new_elements = 0
			old_elements = 0

			for elements, groups in primitive_elements.items():
				if new_pg in groups:
					new_elements = elements
				if old_pg in groups:
					old_elements = elements
			ferro_axes = (int(new_elements) / int(old_elements)) / 2

			os.remove('isotropyout.txt')

			#endregion

			#Sam's new database output - finished

			if interact_with_database == True:   #allows program to run without messing up the database if this is False
				output_list = (ref_code, old_sgname, new_sgname, aizu_yes, continuous, irrep, ferro_axes, run_comments, name, formula, phase_transition, old_spg, new_spg, centro, overall_rss_C, overall_rss_P, overall_rms, max_rms, max_type, max_oldtype, atom_no, sym_group_no, str(temperature), str(melting), str(author), str(journal_name), journal_vol, page_no, year, doi, note, Z, Zprime, calc_density, a_in, b_in, c_in, al_in, be_in, ga_in, vol_in, a_out, b_out, c_out, al_out, be_out, ga_out, vol_out, float(vol_in)/float(vol_out), analogue, bioactivity, ccdc_number, colour, database, depositiondate, disorder, habit, organic, polymeric, powder, polymorph, pressure, r_factor, radiation, remarks, str(solvent), str(source), str(synonyms))
				con = sqlite3.connect('results_database.sqlite')  #connects to database
				cur = con.cursor()
				if new_table == True: #creates new table separate from main 'results' table
					cur.execute("INSERT INTO {tab} VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)".format(tab = new_table_name), output_list)
					con.commit()
					con.close()
				else: #working within the main 'results' table
					cur.execute("SELECT refcode FROM results WHERE refcode = (?)", [ref_code])
					dbrefcode = cur.fetchall()
					dbrefcode = str(dbrefcode)
					dbrefcode = dbrefcode[3:-4]
					if dbrefcode == ref_code: #above code extracts refcode and this checks if it already exists in the database
						new_comment = str(run_version) + ': ' + str(run_comments)
						cur.execute("UPDATE results SET SearchHistory = SearchHistory || char(10) || ? WHERE refcode = ?", (new_comment, ref_code)) #adds run note to final column
						cur.execute("UPDATE results SET Overall_rss_C = ? WHERE refcode = ?", (overall_rss_C, ref_code))
						cur.execute("UPDATE results SET Overall_rss_P = ? WHERE refcode = ?", (overall_rss_P, ref_code))
						entries_updated += 1
						con.commit()
						con.close()
					else:
						cur.execute("INSERT INTO results VALUES (NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", output_list) #adds new entry to 'results' table
						new_comment = run_version + ': ' + str(run_comments)
						cur.execute("UPDATE results SET SearchHistory = ? WHERE refcode = ?",(new_comment, ref_code)) #adds new comment for first time
						entries_added += 1
						con.commit()
						con.close()

			if compare_cifs:
				SUPERIMPOSE.main(ref_code)

		else: #old_sym_no == new_sym_no:  # to save space we delete all the files for irrelevant entries
			if delete_unchanged_cif_inp_out:
				print('No new symmetries found, deleting unnecessary files')
				os.remove(cif)
				os.remove(inp)
				os.remove(out)
				os.remove(findsym_cif)
		print('Number of changed symmetries found:',changed)

	with open('good_ferros.gcd') as infile: #with open('known_ferros.gcd') as infile:
		gcd_list = infile.readlines()
		gcd_list = [x.strip() for x in gcd_list]
	#known_list = ['GUMMUW', 'TEJKUO', 'JAYPUU01', 'UBOTIQ01', 'ZIYYOX', 'ZIYYUD', 'ZIYZAK', 'TTFCAN03', 'TEDAPC01', 'WOLYUR01', 'WOLYUR04', 'MAMPUM02', 'CBUDCX', 'PROLON', 'TAPZIT', 'FAWSAZ01', 'FINTAZ', 'FINTIH', 'FIGNER', 'BENZIL03', 'THIOUR14', 'FOXNUB07', 'JOJYOY', 'LEQWAF', 'YOCVAP', 'YOCVAP01', 'BILNES', 'BILNOC', 'BILPOE01', 'BILNUI', 'HOHMOG01', 'JOJZAL', 'SOWLEW', 'BOBVIY05', 'MIGQOJ02', 'YASKIP', 'BILHUC', 'BILJOY', 'KUWHOA01', 'PIKDOC', 'KASXUY', 'BIKREV', 'FUGMIF', 'ABIQOU', 'ABIRAH', 'AMBACO07', 'ANPHPR01', 'BOLDIP11', 'CYHDAC', 'FOSGID', 'GLYCIN', 'QODGUM04', 'TMPPIO03', 'JIDBIK', 'REZBOP01', 'RIBLET', 'SIYWUT', 'TIVFEL', 'TIVDAF', 'WIVDUA05', 'AMHCLA01', 'DANFAC', 'DGLYCN01', 'DOQVEN01', 'DOQVEN02', 'FEMNEQ01', 'LABSEM03', 'LABSEM04', 'NUKBUS', 'RIDFOA01', 'RIJVIS01', 'RIJVIS02', 'SOMDOO', 'SOMDUU', 'TAXBUQ', 'TGLYBE20', 'TGLYSU', 'XERBAX', 'ZOKKER', 'KOWYEA02', 'DAHTUB', 'HOHMOG05', 'BILGAH', 'BILHAI', 'BILJIS', 'BILKEP', 'BILLEQ', 'BILLOA', 'BILLUG', 'BILMER', 'HORFAV', 'JAWQUT01', 'MEGJEP', 'CASRPP01', 'CATXED01', 'CEKLIN02', 'COHJEP', 'COWBOG', 'DADTOS03', 'DIYJUQ04', 'EGAVAK', 'EHAKII', 'HOJLAU', 'JIKJIW', 'NIXCOP', 'NERDIA', 'NONYUK02', 'SIWKEP', 'OROWAV', 'PIWGEI', 'PIWGIM', 'QOSPAQ', 'QQQBVJ04', 'QQQBVJ31', 'REZXOJ01', 'SARCAC01', 'SATART05', 'SAXYOI', 'TAPBOB06', 'TEDFUF', 'TIXMOD01', 'TMACZN06', 'XUBBUS', 'ZZZUYE04', 'ZZZUYE06', 'UZIWOQ', 'FUWWOL02', 'NIWJUA02', 'WULSAY', 'WULSAY02', 'TOZFAP01', 'TOZFAP02', 'YEYZIN03', 'ELUQIM', 'KNATAR05', 'GLYAGN06', 'GUALSU10', 'GUCRSU', 'ECIDUP', 'WERZUR', 'QODZER', 'OPAJUL01', 'VEMYAN', 'XUQQAB', 'KEFWEX02', 'FIVDAQ', 'FIVDEU', 'EWEBIR', 'EWEBOX', 'CAJXIV', 'AQUPEI01', 'QQQGOP01']
	#['GUMMUW', 'GUHHAT', 'TEJKUO', 'JAYPUU01', 'UBOTIQ01', 'ZIYYOX', 'ZIYYUD', 'ZIYZAK', 'TTFCAN03', 'TEDAPC01', 'WOLYUR01', 'WOLYUR04', 'MAMPUM02', 'CBUDCX', 'PROLON', 'TAPZIT', 'FAWSAZ01', 'FINTAZ', 'FINTIH', 'FIGNER', 'BENZIL03', 'THIOUR14', 'FOXNUB07', 'JOJYOY', 'LEQWAF', 'YOCVAP', 'YOCVAP01', 'BILNES', 'BILNOC', 'BILPOE01', 'BILNUI', 'HOHMOG01', 'JOJZAL', 'SOWLEW', 'BOBVIY05', 'MIGQOJ02', 'YASKIP', 'BILHUC', 'BILJOY', 'KUWHOA01',	 'PIKDOC', 'KASXUY', 'BIKREV', 'FUGMIF', 'ABIQOU', 'ABIRAH', 'AMBACO07', 'ANPHPR01', 'BOLDIP11', 'CYHDAC', 'FOSGID', 'GLYCIN', 'QODGUM04', 'TMPPIO03', 'JIDBIK', 'REZBOP01', 'RIBLET', 'SIYWUT', 'TIVFEL', 'TIVDAF', 'WIVDUA05', 'AMHCLA01', 'DANFAC', 'DGLYCN01', 'DOQVEN01', 'DOQVEN02', 'FEMNEQ01', 'LABSEM03', 'LABSEM04', 'NUKBUS', 'RIDFOA01', 'RIJVIS01', 'RIJVIS02', 'SOMDOO', 'SOMDUU', 'TAXBUQ', 'TGLYBE20', 'TGLYSU', 'XERBAX', 'ZOKKER', 'KOWYEA02', 'DAHTUB', 'HOHMOG05', 'BILGAH', 'BILHAI', 'BILJIS', 'BILKEP', 'BILLEQ', 'BILLOA', 'BILLUG', 'BILMER', 'HORFAV', 'JAWQUT01', 'MEGJEP', 'CASRPP01', 'CATXED01', 'CEKLIN02', 'COHJEP', 'COWBOG', 'DADTOS03', 'DIYJUQ04', 'EGAVAK', 'EHAKII', 'HOJLAU', 'JIKJIW', 'NIXCOP', 'NERDIA',               'SIWKEP', 'OROWAV', 'PIWGEI', 'PIWGIM', 'QOSPAQ', 'QQQBVJ04', 'QQQBVJ31', 'REZXOJ01', 'SARCAC01', 'SATART05', 'SAXYOI', 'TAPBOB06', 'TEDFUF', 'TIXMOD01', 'TMACZN06', 'XUBBUS', 'ZZZUYE04', 'ZZZUYE06', 'UZIWOQ', 'FUWWOL02', 'NIWJUA02', 'WULSAY', 'WULSAY02', 'TOZFAP01', 'TOZFAP02', 'YEYZIN03', 'ELUQIM', 'KNATAR05', 'GLYAGN06', 'GUALSU10', 'GUCRSU', 'ECIDUP', 'WERZUR', 'QODZER', 'OPAJUL01', 'VEMYAN', 'XUQQAB', 'KEFWEX02', 'FIVDAQ', 'FIVDEU', 'EWEBIR', 'EWEBOX', 'CAJXIV', 'AQUPEI01', 'QQQGOP01']
	known_list = gcd_list

	FOM = 0  # figure of merit
	known_hits = []
	for entry in known_list:  # tests how many we got of known ferroelectrics
		if entry in hitlist:
			FOM += 1
			known_hits.append(entry)

	with open('known_found.gcd', 'w') as outfile1:
		for entry in known_list:
			if entry in hitlist:
				outfile1.write(entry + '\n')


	with open('known_notfound.gcd', 'w') as outfile2:
		for entry in known_list:
			if entry not in hitlist:
				outfile2.write(entry + '\n')

	print('Done attempting to process all cif files in list, results in ' + results)  # some final info for the results
	print('Tolerances: latticeTolerance = {} atomicPositionTolerance = {} atomicPositionMaxTolerance = {}'.format(lat_tol,atom_tol,atom_rel_tol))
	print('Number attempted: ' + str(attempted))
	print('Number failed (findsym crashes): ' + str(failed))
	print('Number changed: ' + str(changed))
	print('Number changed which go centrosymmetric: ' + str(changed_centro))
	print('Of {} known FEs, found {}'.format(len(known_list), FOM))
	print('These were:',str(known_hits))
	print('Number of entries added to the database:', str(entries_added)) #added by Sam
	print('Number of entries updated in the database:', str(entries_updated)) #added by Sam

	with open(results, 'a') as outfile:
		outfile.write('Done attempting to process all cif files in list \n'
						'Number attempted: ' + str(attempted) + '\n'
						'Number failed (findsym crashes): ' + str(failed) + '\n'
						'Number changed: ' + str(changed) + '\n' +
						'Of {} known values, found {}\n'.format(len(known_list), FOM) + 'These were:'+str(known_hits)+'\n')


if __name__ == '__main__':
	top_processed = True
	topology_relabel = True
	output_shifts = True
	main()
