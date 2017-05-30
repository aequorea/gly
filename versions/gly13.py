#!/usr/bin/python
# gly13.py -- protein glycosylator
# 2016.10.31 -- John Saeger

input_file = "1yfpA.pdb"				# PDB input file
exposed_file = "1yfpA.exposed_atm.pdb"	# PDB file with solvent exposed atoms

def read_pdb(fn, chain):				# read selected chain in PDB file
	pdb = []							# clean it up a little
	for line in open(fn, 'r'):
		l = line.strip().split()
		if l[0] != 'ATOM' or l[4] != chain: continue
		pdb.append([l[0], l[1], l[2], l[3], l[4], l[5]])
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
	exp = []
	for l in pdb:
		if l[2] not in ['N', 'C', 'O', 'CA', 'CB']:	# more than backbone
			if int(l[5]) not in exp: exp.append(int(l[5]))
	pairs = []							# get solvent exposed pairs
	for s in exp:
		if s+2 in exp: pairs.append(s)
	return pairs

def get_glyco_sites(seq, allowed):
	hydro = { 'A':18, 'R':-45, 'N':-35, 'D':-35, 'C':25, 'E':-35,
		'Q':-35, 'G':-4, 'H':-32, 'I':45, 'L':38, 'K':-39, 'M':19,
		'F':28, 'P':-16, 'S':-8, 'T':-7, 'W':-9, 'Y':-13, 'V':42,
		'Z':0 }

	scores = []							# score potential glycosylation sites ...
	for n in range(2, len(seq)-3):		# ... based on hydropathy
		score = hydro[seq[n-2]]			# hydrophobic
		score += hydro[seq[n-1]]		# hydrophobic
		score -= hydro[seq[n]]			# hydrophilic
		score += hydro[seq[n+1]]		# hydrophobic
		score -= hydro[seq[n+2]]		# hydrophilic
		score += hydro[seq[n+3]]		# hydrophobic
		if 'P' in seq[n-2:n+4]: score = -999	# exclude prolines
		if 'Z' in seq[n-2:n+4]: score = -999	# exclude skipped or unknown residues
		scores.append([n+1, int(score)]) # sequence starts at 0, PDB starts at 1
		n += 1
	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)
	glyco = []								# narrow down the choices
	for sc in sorted_scores:
		if sc[1] < 0: break					# score too low
		if sc[0] not in allowed: continue
		for s in glyco:
			if abs(s[0]-sc[0]) < 3: break	# too close to an existing site
		else:
			delta_h = int(hydro['N']+hydro['T']-hydro[seq[sc[0]-1]]-hydro[seq[sc[0]+1]])
			glyco.append([sc[0], sc[1], seq[sc[0]-1:sc[0]+2], delta_h])
	return glyco

# main

for chain in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
	try:
		sequence = get_sequence(input_file, chain)
		if len(sequence) == 0: continue
		print("\nChain {} sequence has {} residues:\n{}\n".format(chain, len(sequence), sequence))
	except:
		print("Couldn't read input file: {}\n".format(input_file))
		quit()

	try:
		pairs = get_exposed(exposed_file, chain)		# get solvent exposed residue pairs
		print("Chain {} has {} exposed pairs with first residue at:\n{}\n".format(chain, len(pairs), pairs))
	except:
		print("Couldn't read exposed file: {}".format(exposed_file))
		quit()

	glyco_sites = get_glyco_sites(sequence, pairs)
	print("Chain {} has {} potential glycosylation sites:".format(chain, len(glyco_sites)))
	print("[site, score, seq, delta hydropathy]")
	for gs in glyco_sites: print gs
