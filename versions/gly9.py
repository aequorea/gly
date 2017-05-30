#!/usr/bin/python
# gly9.py -- protein glycosylator
# John Saeger (2016.10.25)

prefix = "1yfpA"						# file prefix of molecule to work on

input_file = prefix+".pdb"				# PDB input file
exclude_file = prefix+".exclude"		# file with excluded residues
exposed_file = prefix+".exposed_atm.pdb" # PDB file with solvent exposed atoms

glyco_file = prefix+".glyco"			# gets list of glycosylation sites
score_file = prefix+".score"			# gets list of glycosylation scores


def read_fasta(fn):						# get strings -- ignore comment lines
	return ''.join(["" if l[0] in ";>" else l.strip() for l in open(fn, 'r')])
		
def read_pdb(fn):						# get strings
	return [l.strip().split() for l in open(fn, 'r')]
	
def read_string_list(fn):				# get strings -- convert to ints
	return [int(l.strip()) for l in open(fn, 'r')]
		
def write_int_list(fn, glyco):			# write ints as strings
	with open(fn, "w") as fh: fh.write('\n'.join(str(x) for x in glyco)+'\n')
	return
	
def	get_seq(pdb):						# get fasta style sequence from PDB file
										# pad skips and unknown with 'Z'
										
	res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
		
	seq = ""
	resnum = 0
	for l in pdb:
		if l[0] != 'ATOM': continue
		if l[4] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ": # see if we have a chain name
			newresnum = int(l[5])
		else:
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
	
def get_exposed(pdb):					# get list of residues with ...
										# ... solvent exposed sidechains
	first = []
	exp = []
	for l in pdb:
		if l[0] != 'ATOM': continue
		if l[4] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ": # see if we have a chain name
			newresnum = int(l[5])
		else:
			newresnum = int(l[4])
		if l[2] not in ['N', 'C', 'O', 'CA', 'CB']:	# beyond backbone
			if newresnum not in first:
				first.append(newresnum)
				
	for s in first:
		if s+2 in first:				# exposed residue is paired
			exp.append(s)
			
	return(first, exp)

def get_glyco_sites(sequence, exclude):
	
	hydropathy = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2,
		'Z':0 }

	scores = []								# score potential glycosylation sites
	for n in range(2, len(sequence)-3):
		score = hydropathy[sequence[n-2]]*10		# hydrophobic
		score += hydropathy[sequence[n-1]]*10		# hydrophobic
		score -= hydropathy[sequence[n]]*10			# hydrophilic
		score += hydropathy[sequence[n+1]]*10		# hydrophobic
		score -= hydropathy[sequence[n+2]]*10		# hydrophilic
		score += hydropathy[sequence[n+3]]*10		# hydrophobic
		if 'P' in sequence[n-2:n+4]: score = -999	# exclude prolines
		if 'Z' in sequence[n-2:n+4]: score = -999	# exclude skipped or unknown residues
		scores.append([n, int(score)])
		n += 1

	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)

	glyco = []								# narrow down the choices
	for sc in sorted_scores:
		site, score = sc[0], sc[1]
		if score < 0: break
		if site in exclude: continue
		for s in glyco:
			if abs(s[0]-site) < 3: break	# too close to an existing site
		else:
			glyco.append([site, score, sequence[site-1]])

	return glyco
	

# main

try:
	pdb = read_pdb(input_file)
except:
	print("\nCouldn't read PDB file: {}".format(pdb_file))
	quit()

print("PDB file has {} lines.".format(len(pdb)))

sequence = get_seq(pdb)
print("Sequence has {} residues:".format(len(sequence)))
print(sequence)

try:
	exclude = read_string_list(exclude_file)
except:

	print("\nCouldn't read exclude file: {}".format(exclude_file))
	print("Assuming no manual exclusions.")
	exclude = []

print("\nManual exclusions:")	
print(exclude)
	
try:
	exposed = read_pdb(exposed_file)
except:
	print("\nCouldn't read exposed file: {}".format(exposed_file))
	print("Assuming full solvent exposure.\n")
else:
	first, exp = get_exposed(exposed)
	print("\nThere are {} residues with exposed side chains:".format(len(first)))
	print(first)
	print("\nThere are {} exposed residues with a second exposed residue at n+2:".format(len(exp)))
	print(exp)
	for r in range(1,len(sequence)+1):
		if r in exclude: continue
		if r not in exp:
			exclude.append(r)
			
print("\nManual and solvent exclusions:")	
print(exclude)
sites, scores, res = zip(*get_glyco_sites(sequence, exclude))

write_int_list(glyco_file, sites)
write_int_list(score_file, scores)

print("\nThere are {} potential glycosylation sites (site, score):".format(len(sites)))
print(zip(sites, scores, res))

print(sites)
print
print(scores)

