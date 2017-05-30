#!/usr/bin/python
# gly5.py -- protein glycosylator
# John Saeger

infile = "lasp.fas"			# fasta input file
num_sites = 16				# desired number of glycosylation sites
exclusions = []				# list of sites to disallow

def read_fasta(fn):
	fasta_str = ""
	seq = ""
	with open(fn) as fh:
		for line in fh:		# sort out the comment lines
			if line[0] == '>' or line[0] == ';':
				fasta_str += line
			else:
				seq += line.rstrip()
	return fasta_str, seq	# return comment lines and sequence
			
def score_site(seq, n):
	
	hydropathy = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
		'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
		'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2 }

	if 'P' in seq[n-2:n+4]: 
		return -999

	score = 0
	score += hydropathy[seq[n-2]]*10	# plus  - hydrophobic
	score += hydropathy[seq[n-1]]*10	# plus  - hydrophobic
	score -= hydropathy[seq[n]]*10		# minus - hydrophilic
	score += hydropathy[seq[n+1]]*10	# plus  - hydrophobic
	score -= hydropathy[seq[n+2]]*10	# minus - hydrophilic
	score += hydropathy[seq[n+3]]*10	# plus  - hydrophobic	
	return int(score)
	
def get_glyco_sites(infile, num_sites, exclusions):
	fasta_str, sequence = read_fasta(infile)
	lower_limit = 2						# where to look for glycosylation sites
	upper_limit = len(sequence) - 3
	
	n = lower_limit
	scores = []
	while n < upper_limit:				# score possible glycosylation sites
		site_score = score_site(sequence, n)
		scores.append([n, site_score])
		n += 1

	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)

	good_sites = 0
	top_sites = []
	glyco = []							
	glyco_scores = []					
	too_close = []
	too_close_scores = []
	for sc in sorted_scores:
		top = sc[0]
		score = sc[1]
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

	return glyco, glyco_scores, too_close, too_close_scores
	
# main
	
gly, gly_scores, tc, tc_scores = get_glyco_sites(infile, num_sites, exclusions)
print(gly)
print(gly_scores)
print(tc)
print(tc_scores)
