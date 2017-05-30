#!/usr/bin/python
# gly8.py -- protein glycosylator
# John Saeger (2016.10.25)


num_sites = 10							# number of sites to find -- (-1 = all)
prefix = "1yfpA"						# file prefix of molecule to work on

fasta_file = prefix+".fasta"			# fasta input file
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
										# pad skips with 'Z'
										
	res_letters = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',
		'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
		'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
	
	seq = ""
	resnum = 0
	for l in pdb:
		if l[0] != 'ATOM': continue
		if int(l[5]) != resnum:
			if int(l[5]) == resnum+1:
				resnum += 1
				resname = l[3]
				seq += res_letters[resname]
			else:
				resnum += 1
				seq += 'Z'
	return seq
	
def get_exposed(pdb):
	exp = []
	for l in pdb:
		if l[0] != 'ATOM': continue
		if l[2] not in ['N', 'CA', 'C', 'O']:
			if int(l[5]) not in exp:
				exp.append(int(l[5]))
		if l[3] == 'GLY' and l[2] == 'CA':
			exp.append(int(l[5]))
	return(exp)

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
			glyco.append([site, score])

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
	print("Assuming no manual exclusions.\n")
	exclude = []

print("Manual exclusions:")	
print(exclude)
	
try:
	exposed = read_pdb(exposed_file)
except:
	print("\nCouldn't read exposed file: {}".format(exposed_file))
	print("Assuming full solvent exposure.\n")
else:
	exp = get_exposed(exposed)
	print("\nThere are {} exposed residues:".format(len(exp)))
	print(exp)
	for r in range(1,len(sequence)+1):
		if r in exclude: continue
		if r not in exp:
			exclude.append(r)
			
print("\nManual and solvent exclusions:")	
print(exclude)
sites, scores = zip(*get_glyco_sites(sequence, exclude))

write_int_list(glyco_file, sites[:num_sites])
write_int_list(score_file, scores[:num_sites])

print("\nGlycosylation (site, score):")
print(zip(sites, scores)[:num_sites])
