#!/usr/bin/python
# gly18.py -- protein glycosylator
# 2016.11.03 -- John Saeger

file_prefix = "1yfp"

input_file = file_prefix + ".pdb"				# PDB input file
exposed_file = file_prefix + "_exposed_atm.pdb"	# PDB file with solvent exposed atoms
sheet_file = file_prefix + "_sheet.pdb"			# PDB file with beta sheet residues
helix_file = file_prefix + "_helix.pdb"			# PDB file with alpha helix residues

def read_pdb(fn, chain):				# read selected chain in PDB file
	try:
		pdb = []							# clean it up a little
		for line in open(fn, 'r'):
			l = line.strip().split()
			if l[0] != 'ATOM' or l[4] != chain: continue
			pdb.append([l[0], l[1], l[2], l[3], l[4], l[5]])
	except:
		print("Couldn't read {}".format(fn))
		quit()
	return pdb

def	get_sequence(fn, chain):			# get sequence of selected chain from PDB file
	res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

	pdb = read_pdb(fn, chain)
	seq = ""
	resnum = 0
	for l in pdb:
		if int(l[5]) != resnum:
			if int(l[5]) == resnum+1:	# see if there is a resnum skip
				resnum += 1
				try:
					seq += res_letters[l[3]]
				except:
					seq += 'Z'			# set to 'Z' if residue unknown
			else:
				resnum += 1
				seq += 'Z'				# set to 'Z' if there was a skip
	return seq

def get_exposed(fn, chain):				# get solvent exposed residues
	pdb = read_pdb(fn, chain)
	singles = []
	for l in pdb:
		if l[2] not in ['N', 'C', 'O', 'CA', 'CB']:	# more than backbone
			if int(l[5]) not in singles: singles.append(int(l[5]))
	pairs = []							# get solvent exposed pairs
	for s in singles:
		if s+2 in singles: pairs.append(s)
	return pairs
	
def get_residues(fn, chain):
	pdb = read_pdb(fn, chain)
	residues = []
	for l in pdb:
		if int(l[5]) not in residues:
			residues.append(int(l[5]))
	return residues
	
def get_random_coil(input_file, chain, sheet, helix):
	all = get_residues(input_file, chain)
	coil = []
	for res in all:
		if res not in sheet and res not in helix:
			coil.append(res)
	return coil

def sheet_score(hydro, seq, n):
	score = 0
	score += hydro[seq[n-2]]		# hydrophobic
	score += hydro[seq[n-1]]		# hydrophobic
	score -= hydro[seq[n]]			# hydrophilic
	score += hydro[seq[n+1]]		# hydrophobic
	score -= hydro[seq[n+2]]		# hydrophilic
	score += hydro[seq[n+3]]		# hydrophobic
	return score
	
def helix_score(hydro, seq, n):
	score = 0
#	score -= hydro[seq[n-2]]		# hydrophilic
#	score -= hydro[seq[n-1]]		# hydrophilic
	score -= hydro[seq[n]]			# hydrophilic
	score -= hydro[seq[n+1]]		# hydrophilic
	score -= hydro[seq[n+2]]		# hydrophilic
#	score -= hydro[seq[n+3]]		# hydrophilic
	
	return score

def get_glyco_sites(score_func, seq, allowed1, allowed2):
	hydro = { 'A':18, 'R':-45, 'N':-35, 'D':-35, 'C':25, 'E':-35,
		'Q':-35, 'G':-4, 'H':-32, 'I':45, 'L':38, 'K':-39, 'M':19,
		'F':28, 'P':-16, 'S':-8, 'T':-7, 'W':-9, 'Y':-13, 'V':42,
		'Z':0 }

	scores = []							# score potential glycosylation sites ...
	for n in range(2, len(seq)-3):		# ... based on hydropathy
		score = score_func(hydro, seq, n)
		if seq[n+1] == 'N': score = -999		# exclude N in the middle
		if 'P' in seq[n-2:n+4]: score = -999	# exclude prolines
		if 'Z' in seq[n-2:n+4]: score = -999	# exclude skipped or unknown residues
		scores.append([n+1, int(score)]) # sequence starts at 0, PDB starts at 1
		n += 1
	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)
	glyco = []								# narrow down the choices
	glyco_short = []
	for sc in sorted_scores:
		if sc[1] <= 0: break					# score too low
		if sc[0] not in allowed1 or sc[0] not in allowed2: continue
		for s in glyco:
			if abs(s[0]-sc[0]) < 3: break	# too close to an existing site
		else:
			delta_h = int(hydro['N']+hydro['T']-hydro[seq[sc[0]-1]]-hydro[seq[sc[0]+1]])
			glyco.append([sc[0], sc[1], seq[sc[0]-1:sc[0]+2], delta_h])
			glyco_short.append(sc[0])
	return glyco, glyco_short
	
def show_sites(common, score_func, name, seq, allowed1, allowed2):
	sites, short = get_glyco_sites(score_func, seq, allowed1, allowed2)
	print("\nChain {} has {} potential {} glycosylation sites:".format(chain, len(sites), name))
	print("[site, score, seq, delta hydropathy]")
	for gs in sites: print gs

	if len(common) == 0:
		common.update(short)
	else:
		common = common & set(short)
	return common

# main

sheet_common = set()
helix_common = set()
coil_common = set()
for chain in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
	sequence = get_sequence(input_file, chain)
	if len(sequence) == 0: continue
	print("\nChain {} sequence has {} residues:\n{}".format(chain, len(sequence), sequence))
	pairs = get_exposed(exposed_file, chain)
	sheet = get_residues(sheet_file, chain)
	helix = get_residues(helix_file, chain)
	coil = get_random_coil(input_file, chain, sheet, helix)
	
	sheet_common = show_sites(sheet_common, sheet_score, "beta sheet", sequence, pairs, sheet)
	helix_common = show_sites(helix_common, helix_score, "alpha helix", sequence, pairs, helix)
	coil_common = show_sites(coil_common, helix_score, "random coil", sequence, pairs, coil)
	
print("\nThere are {} common potential beta sheet glycosylation sites:".format(len(sheet_common)))
print(list(sheet_common))

print("\nThere are {} common potential alpha helix glycosylation sites:".format(len(helix_common)))
print(list(helix_common))

print("\nThere are {} common potential random coil glycosylation sites:".format(len(coil_common)))
print(list(coil_common))
