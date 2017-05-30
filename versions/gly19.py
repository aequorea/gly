#!/usr/bin/python

# gly19.py -- protein glycosylator
# 2016.10.05 -- first version -- John Saeger

# gly19 -- 2017.5.27

# This version complains about incomplete sidechains. You must fix them up.
# I use pymol's mutate wizard to mutate a residue to itself.

# This version has a slightly more sophisticated solvent exposure test.
# You can tune it with the solv_thresh variable (below).

# This version doesn't check if glycosylation sites are too close to 
# each other, so it prints out more sites. Use judgement when placing
# glycans close together. You probably won't want to try all sites at
# the same time. It also doesn't cut off based on the bioinformatics
# scoring algorithm. We show everything except things like sites with
# prolines in the immediate neighborhood. After all, I don't know for sure
# if the bioinformatics based scoring is actually telling me anything.

file_prefix = "1yfpA"

input_file = file_prefix + ".pdb"				# PDB input file
exposed_file = file_prefix + "_exposed_atm.pdb"	# PDB file with solvent exposed atoms
sheet_file = file_prefix + "_sheet.pdb"			# PDB file with beta sheet residues
helix_file = file_prefix + "_helix.pdb"			# PDB file with alpha helix residues

solv_thresh = 0.33	# this fraction of sidechain atoms exposed to declare residue solvent exposed

res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

res_atoms = {'ALA':6, 'ARG':12, 'ASN':9, 'ASP':9, 'CYS':7, 'GLU':10,
		'GLN':10, 'GLY':5, 'HIS':11, 'ILE':9, 'LEU':9, 'LYS':10, 'MET':9,
		'PHE':12, 'PRO':8, 'SER':7, 'THR':8, 'TRP':15, 'TYR':13, 'VAL':8}
		
hydro = { 'A':18, 'R':-45, 'N':-35, 'D':-35, 'C':25, 'E':-35,
		'Q':-35, 'G':-4, 'H':-32, 'I':45, 'L':38, 'K':-39, 'M':19,
		'F':28, 'P':-16, 'S':-8, 'T':-7, 'W':-9, 'Y':-13, 'V':42,
		'Z':0 }

def read_pdb(fn, chain):				# read selected chain in PDB file
	try:
		pdb = []						# clean it up a little
		for line in open(fn, 'r'):
			l = line.strip().split()
			if l[0] != 'ATOM' or l[4] != chain: continue
			pdb.append([l[0], l[1], l[2], l[3], l[4], l[5]])
	except:
		print("Couldn't read {}".format(fn))
		quit()
	return pdb

def	get_sequence(fn, chain):			# get sequence of selected chain from PDB file
	pdb = read_pdb(fn, chain)
	seq = ""
	resnum = 0
	atomnum = 0
	for l in pdb:
		if resnum == 0 and atomnum == 0:
			l=pdb[0]
			resname = l[3]
		atomnum += 1
		if int(l[5]) != resnum:
					
			oldresnum = resnum
			oldresname = resname
			oldatomnum = atomnum
								
			resname = l[3]
			atomnum = 0
			if int(l[5]) != resnum+1:	# see if there is a resnum skip
				while resnum < int(l[5])-1:
					resnum += 1
					seq += 'Z'			# set to 'Z' if there was a skip
			resnum += 1
			try:
				seq += res_letters[l[3]]
				if res_atoms[oldresname]-1 != oldatomnum and oldresnum != 0:
					print("Warning! residue {} {} needs {} atoms but only has {}").format(oldresnum, oldresname, res_atoms[oldresname]-1, oldatomnum)
			except:
				seq += 'Z'				# set to 'Z' if residue unknown
	return seq

def	get_singles(fn, chain):				# get solvent exposed atoms from PDB file
	pdb = read_pdb(fn, chain)
	resnum = 0
	atoms = 0							# we'll be counting sidechain atoms
	singles = []

	for l in pdb:
		if resnum == 0 and atoms == 0:
			l=pdb[0]
			resname = l[3]
			
		if l[2] not in ['N', 'C', 'O'] and int(l[5]) == resnum:		# had CA CB
			atoms += 1
		if int(l[5]) != resnum:
					
			oldresnum = resnum
			oldresname = resname
			oldatoms = atoms
								
			resname = l[3]
			resnum = int(l[5])
			
			if l[2] not in ['N', 'C', 'O']:
				atoms = 1
			else:
				atoms = 0
				
			if oldresname in res_atoms:
				if oldatoms > solv_thresh * (res_atoms[oldresname]-4) and oldresnum != 0:
					singles.append(oldresnum)
	return singles

def get_exposed(fn, chain):				# get solvent exposed residues
	pdb = read_pdb(fn, chain)

	singles = get_singles(fn, chain)

	pairs = []							# get solvent exposed pairs
	for s in singles:
		if s+2 in singles: pairs.append(s)
	return singles, pairs
	
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
	score += hydro[seq[n-2]]			# hydrophobic
	score += hydro[seq[n-1]]			# hydrophobic
	score -= hydro[seq[n]]				# hydrophilic
	score += hydro[seq[n+1]]			# hydrophobic
	score -= hydro[seq[n+2]]			# hydrophilic
	score += hydro[seq[n+3]]			# hydrophobic
	return score
	
def helix_score(hydro, seq, n):
	score = 0
#	score -= hydro[seq[n-2]]			# hydrophilic
#	score -= hydro[seq[n-1]]			# hydrophilic
	score -= hydro[seq[n]]				# hydrophilic
	score -= hydro[seq[n+1]]			# hydrophilic
	score -= hydro[seq[n+2]]			# hydrophilic
#	score -= hydro[seq[n+3]]			# hydrophilic
	
	return score

def get_glyco_sites(score_func, seq, allowed1, allowed2):
	scores = []									# score potential glycosylation sites ...
	for n in range(2, len(seq)-3):				# ... based on hydropathy
		score = score_func(hydro, seq, n)
#		if seq[n+1] == 'N': score = -999		# exclude N in the middle
		if 'P' in seq[n-2:n+4]: score = -999	# exclude prolines
		if 'Z' in seq[n-2:n+4]: score = -999	# exclude skipped or unknown residues
		scores.append([n+1, int(score)]) 		# sequence starts at 0, PDB starts at 1
		n += 1
#	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)
	glyco = []									# narrow down the choices
	glyco_short = []
	for sc in scores:
		if sc[1] <= -999: continue				# score too low
		if sc[0] not in allowed1 or sc[0] not in allowed2: continue
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
		
	scommon = sorted(common)
	return scommon

# main

sheet_common = set()
helix_common = set()
coil_common = set()
total_common = set()
for chain in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
	sequence = get_sequence(input_file, chain)
	if len(sequence) == 0: continue
	print("\nChain {} sequence has {} residues:\n{}".format(chain, len(sequence), sequence))
	
	singles, pairs = get_exposed(exposed_file, chain)
	print("\nChain {} has {} solvent accessible residues:\n{}".format(chain, len(singles), singles))
	print("\nChain {} has {} solvent accessible pairs:\n{}".format(chain, len(pairs), pairs))
	
	sheet = get_residues(sheet_file, chain)
	print("\nChain {} has {} beta sheet residues:\n{}".format(chain, len(sheet), sheet))
	helix = get_residues(helix_file, chain)
	print("\nChain {} has {} alpha helix residues:\n{}".format(chain, len(helix), helix))
	coil = get_random_coil(input_file, chain, sheet, helix)
	print("\nChain {} has {} random coil residues:\n{}".format(chain, len(coil), coil))
	
	sheet_common = show_sites(sheet_common, sheet_score, "beta sheet", sequence, pairs, sheet)
	helix_common = show_sites(helix_common, helix_score, "alpha helix", sequence, pairs, helix)
	coil_common = show_sites(coil_common, helix_score, "random coil", sequence, pairs, coil)
	
print("\nThere are {} common potential beta sheet glycosylation sites:".format(len(sheet_common)))
print(list(sheet_common))

print("\nThere are {} common potential alpha helix glycosylation sites:".format(len(helix_common)))
print(list(helix_common))

print("\nThere are {} common potential random coil glycosylation sites:".format(len(coil_common)))
print(list(coil_common))

total_common = sheet_common + helix_common + coil_common
stotal_common = sorted(list(total_common))

print("\nThere are {} common potential glycosylation sites:".format(len(stotal_common)))
print(stotal_common)

