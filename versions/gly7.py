#!/usr/bin/python
# gly7.py -- protein glycosylator
# John Saeger (2016.10.25)


num_sites = 10							# number of sites to find -- (-1 = all)
prefix = "eyfp"							# file prefix of molecule to work on

fasta_file = prefix+".fasta"			# fasta input file
exclude_file = prefix+".exclude"		# file with excluded residues
solvent_file = prefix+".solvent"		# file with solvent accessible residues

glyco_file = prefix+".glyco"			# gets list of glycosylation sites
score_file = prefix+".score"			# gets list of glycosylation scores


def read_fasta(fn):						# get strings -- ignore comment lines
	return ''.join(["" if l[0] in ";>" else l.strip() for l in open(fn, 'r')])
		
def read_exclude(fn):					# get strings -- convert to ints
	return [int(l.strip()) for l in open(fn, 'r')]
		
def write_int_list(fn, glyco):			# write ints as strings
	with open(fn, "w") as fh: fh.write('\n'.join(str(x) for x in glyco)+'\n')
	return

def get_glyco_sites(sequence, exclude):
	
	hydropathy = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2 }

	scores = []								# score potential glycosylation sites
	for n in range(2, len(sequence)-3):
		score = hydropathy[sequence[n-2]]*10		# hydrophobic
		score += hydropathy[sequence[n-1]]*10		# hydrophobic
		score -= hydropathy[sequence[n]]*10			# hydrophilic
		score += hydropathy[sequence[n+1]]*10		# hydrophobic
		score -= hydropathy[sequence[n+2]]*10		# hydrophilic
		score += hydropathy[sequence[n+3]]*10		# hydrophobic
		if 'P' in sequence[n-2:n+4]: score = -999	# exclude prolines
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
	sequence = read_fasta(fasta_file)
except:
	print("Couldn't read fasta file: {}\n".format(fasta_file))
	quit()

print("\nSequence:")
print(sequence + '\n')
	
try:
	exclude = read_exclude(exclude_file)
except:
	print("Couldn't read exclude file: {}".format(exclude_file))
	print("Assuming no exclusions.\n")
	exclude = []

print("Exclusions:")	
print(exclude)
	
sites, scores = zip(*get_glyco_sites(sequence, exclude))

write_int_list(glyco_file, sites[:num_sites])
write_int_list(score_file, scores[:num_sites])

print("\nGlycosylation (site, score):")
print(zip(sites, scores)[:num_sites])
