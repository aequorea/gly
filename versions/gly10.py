#!/usr/bin/python
# gly10.py -- protein glycosylator
# John Saeger (2016.10.26)

prefix = "1yfpA"						# file prefix of molecule to work on

input_file = prefix+".pdb"				# PDB input file
exclude_file = prefix+".exclude"		# file with excluded residues
exposed_file = prefix+".exposed_atm.pdb" # PDB file with solvent exposed atoms

def read_string_list(fn):				# get strings -- convert to ints
	return [int(l.strip()) for l in open(fn, 'r')]
		
def read_pdb(fn):						# read PDB file -- clean it up a little
	clean = []
	for line in open(fn, 'r'):
		l = line.strip().split()
		if l[0] != 'ATOM': continue
		if l[4] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ": # see if we have a chain name
			clean.append([l[0], l[1], l[2], l[3], l[5]])
		else:
			clean.append([l[0], l[1], l[2], l[3], l[4]])
	return clean

def	get_sequence(fn):					# get fasta style sequence from PDB file										
	res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
		
	pdb = read_pdb(fn)
	seq = ""
	resnum = 0
	for l in pdb:
		newresnum = int(l[4])
		if newresnum != resnum:
			if newresnum == resnum+1:	# see if there is a resnum skip
				resnum += 1
				resname = l[3]
				try:
					seq += res_letters[resname]
				except:
					seq += 'Z'			# set to 'Z' if residue unknown
			else:
				resnum += 1
				seq += 'Z'				# set to 'Z' if there was a skip
	return seq
	
def get_exposed(fn):					# get solvent exposed residues
	pdb = read_pdb(fn)	
	first = []
	exp = []
	for l in pdb:
		newresnum = int(l[4])
		if l[2] not in ['N', 'C', 'O', 'CA', 'CB']:	# beyond backbone
			if newresnum not in first:
				first.append(newresnum)
	for s in first:
		if s+2 in first:				# exposed residue is paired
			exp.append(s)
	return(first, exp)

def get_glyco_sites(seq, exclude):
	hydro = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2,
		'Z':0 }
		
	scores = []							# score potential glycosylation sites
	for n in range(2, len(seq)-3):
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
		if sc[0] in exclude: continue
		for s in glyco:
			if abs(s[0]-sc[0]) < 3: break	# too close to an existing site
		else:
			delta_h = int(10*(hydro['N']+hydro['T']-hydro[seq[sc[0]-1]]-hydro[seq[sc[0]+1]]))
			glyco.append([sc[0], sc[1], seq[sc[0]-1:sc[0]+2], delta_h])
	return glyco
	

# main

try:
	sequence = get_sequence(input_file)
	print("Sequence has {} residues:".format(len(sequence)))
	print(sequence)
except:
	print("\Couldn't read input file: {}".format(input_file))
	quit()

try:
	exclude = read_string_list(exclude_file)
except:
	print("\nCouldn't read exclude file: {}".format(exclude_file))
	print("Assuming no manual exclusions.")
	exclude = []
print("\nManual exclusions:")	
print(exclude)

try:
	first, pairs = get_exposed(exposed_file)
	print("\nThere are {} residues with exposed side chains:".format(len(first)))
	print(first)
	print("\nThere are {} exposed residue pairs:".format(len(pairs)))
	print(pairs)
	for r in range(1,len(sequence)+1):	# add unpaired residues to exclude list
		if r in exclude: continue
		if r not in pairs:
			exclude.append(r)
except:
	print("\nCouldn't read exposed file: {}".format(exposed_file))
	print("Assuming full solvent exposure.\n")

print("\nManual and solvent exclusions:")	
print(exclude)
glyco_sites = get_glyco_sites(sequence, exclude)
print("\nThere are {} potential glycosylation sites (site, score):".format(len(glyco_sites)))
for gs in glyco_sites: print gs
