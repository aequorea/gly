#!/usr/bin/python
# gly6.py -- protein glycosylator
# John Saeger

infile = "lasp.fas"			# fasta input file
exclusions = []				# list of sites to disallow

def read_fasta(fn):
	seq = ""
	with open(fn) as fh:
		for line in fh:		# ignore comment lines
			if line[0] != '>' and line[0] != ';':
				seq += line.rstrip()
	
	return seq	# return comment lines and sequence
			
def get_glyco_sites(sequence, exclusions):
	
	hydropathy = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2 }

	num_sites = 16					# return this many glycosylation sites

	lower_limit = 2					# where to look for glycosylation sites
	upper_limit = len(sequence) - 3
	
	scores = []
	n = lower_limit
	while n < upper_limit:			# score possible glycosylation sites
		
		score = 0
		score += hydropathy[sequence[n-2]]*10	# plus  - hydrophobic
		score += hydropathy[sequence[n-1]]*10	# plus  - hydrophobic
		score -= hydropathy[sequence[n]]*10		# minus - hydrophilic
		score += hydropathy[sequence[n+1]]*10	# plus  - hydrophobic
		score -= hydropathy[sequence[n+2]]*10	# minus - hydrophilic
		score += hydropathy[sequence[n+3]]*10	# plus  - hydrophobic
		
		if 'P' in sequence[n-2:n+4]: score = -999

		scores.append([n, int(score)])
		n += 1

	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)

	good_sites = 0
	top_sites = []
	glyco = []							
	glyco_scores = []					
	too_close = []
	too_close_scores = []
	for sc in sorted_scores:
		top, score = sc[0], sc[1]
		for site in top_sites:			
			if abs(top-site) < 3:		
				too_close.append(top)	# collect sites too close to a ...
				too_close_scores.append(score)	# ... higher scoring site
		top_sites.append(top)			# collect top sites so far
		if good_sites >= num_sites: break
		if top not in exclusions and top not in too_close:
			good_sites += 1
			glyco.append(top)			# collect glycosylation sites ...
			glyco_scores.append(score)	# ... and their scores

	return glyco, glyco_scores
	
# main

sequence = read_fasta(infile)
gly, gly_scores = get_glyco_sites(sequence, exclusions)

print(gly)
print(gly_scores)
