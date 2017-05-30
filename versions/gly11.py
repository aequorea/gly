#!/usr/bin/python
# gly11.py -- protein glycosylator
# John Saeger (2016.10.26)

input_file = "1yfpA.pdb"				# PDB input file
exposed_file = "1yfpA.exposed_atm.pdb"	# PDB file with solvent exposed atoms

def read_pdb(fn):						# read PDB file -- clean it up a little
	pdb = []
	for line in open(fn, 'r'):
		l = line.strip().split()
		if l[0] != 'ATOM': continue
		if l[4] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ": # see if we have a chain name
			pdb.append([l[0], l[1], l[2], l[3], l[5]])
		else:
			pdb.append([l[0], l[1], l[2], l[3], l[4]])
	return pdb

def	get_sequence(fn):					# get sequence from PDB file
	res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

	pdb = read_pdb(fn)
	seq = ""
	resnum = 0
	for l in pdb:
		if int(l[4]) != resnum:
			if int(l[4]) == resnum+1:	# see if there is a resnum skip
				resnum += 1
				try:
					seq += res_letters[l[3]]
				except:
					seq += 'Z'			# set to 'Z' if residue unknown
			else:
				resnum += 1
				seq += 'Z'				# set to 'Z' if there was a skip
	return seq

def get_exposed(fn):					# get solvent exposed residues
	pdb = read_pdb(fn)
	exp = []
	for l in pdb:
		if l[2] not in ['N', 'C', 'O', 'CA', 'CB']:	# more than backbone
			if int(l[4]) not in exp: exp.append(int(l[4]))
	pairs = []							# get solvent exposed pairs
	for s in exp:
		if s+2 in exp: pairs.append(s)
	return pairs

def get_glyco_sites(seq, allowed):
	hydro = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2,
		'Z':0 }

	scores = []							# score potential glycosylation sites ...
	for n in range(2, len(seq)-3):		# ... based on hydropathy
		score = hydro[seq[n-2]]*10		# hydrophobic
		score += hydro[seq[n-1]]*10		# hydrophobic
		score -= hydro[seq[n]]*10		# hydrophilic
		score += hydro[seq[n+1]]*10		# hydrophobic
		score -= hydro[seq[n+2]]*10		# hydrophilic
		score += hydro[seq[n+3]]*10		# hydrophobic
		if 'P' in seq[n-2:n+4]: score = -999	# exclude prolines
		if 'Z' in seq[n-2:n+4]: score = -999	# exclude skipped or unknown residues
		scores.append([n, int(score)])
		n += 1
	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)
	glyco = []								# narrow down the choices
	for sc in sorted_scores:
		if sc[1] < 0: break					# score too low
		if sc[0] not in allowed: continue
		for s in glyco:
			if abs(s[0]-sc[0]) < 3: break	# too close to an existing site
		else:
			delta_h = int(10*(hydro['N']+hydro['T']-hydro[seq[sc[0]-1]]-hydro[seq[sc[0]+1]]))
			glyco.append([sc[0], sc[1], seq[sc[0]-1:sc[0]+2], delta_h])
	return glyco

# main

try:
	sequence = get_sequence(input_file)
	print("Sequence has {} residues:\n{}\n".format(len(sequence), sequence))
except:
	print("Couldn't read input file: {}\n".format(input_file))
	quit()

try:
	pairs = get_exposed(exposed_file)		# get solvent exposed residue pairs
	print("There are {} exposed pairs with first residue at:\n{}\n".format(len(pairs), pairs))
except:
	print("Couldn't read exposed file: {}".format(exposed_file))
	quit()

glyco_sites = get_glyco_sites(sequence, pairs)
print("There are {} potential glycosylation sites:".format(len(glyco_sites)))
print("[site, score, seq, delta hydropathy]")
for gs in glyco_sites: print gs
