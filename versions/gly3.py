
# gly2.py -- protein glycosylator
# 2016-10-05 -- first version (John Saeger)
# 2016-10-06 -- second version (JS)

# This version attempts to implement what is known about glycosylation
# sites as reported in this study: 

# The Structural Assessment of Glycosylation Sites Database - SAGS -
# An Overall View on N-Glycosylation by Marius D. Surleac, 
# Laurentiu N. Spiridon, Robi Tacutu, Adina L Milac, 
# Stefana M. Petrescu and Andrei-J Petrescu 
# (http://dx.doi.org/10.5772/51690)

# The version 2 scoring algorithm uses the same information but in a 
# more abstract way. We use the pattern from n-2 to n+3 of
# hydrophobic, hydrophobic, hydrophilic, hydrophobic, hydrophilic, 
# hydrophobic from the above study, even though I don't understand why 
# there are two initial hydrophobic residues. We use the hydropathy 
# index from Wikipedia (https://en.wikipedia.org/wiki/Amino_acid)
# to score the pattern. Prolines are excluded and we insist on seeing
# a complete pattern, thus we don't glycosylate near the very ends of 
# the protein. We also don't allow a glycosylation site within 2 or 
# less positions of a higher scoring site. This prevents glycosylation 
# sites from trampling each other.


infile = "lasp.fas"			# fasta input file
outfile = "lasp.gly2"		# second version fasta output
dumpfile = "lasp.dump"		# dump all scores
num_sites = 20				# desired number of glycosylation sites
exclusions = []				# list of sites to disallow

hydropathy = { 'A':1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C':2.5, 'E':-3.5,
	'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I':4.5, 'L':3.8, 'K':-3.9, 'M':1.9,
	'F':2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V':4.2 }
	

#scores = []
#top_sites = []
#top_scores = []
#too_close = []
#gly = ""

def read_fasta(fn):
	fasta_str = ""
	seq = ""
	with open(fn) as fh:
		for line in fh:
			if line[0] == '>' or line[0] == ';':
				fasta_str += line
			else:
				seq += line.rstrip()
	return [fasta_str, seq]
	
def format_score(loc, val, seq, exclusions, too_close):
	ss = seq[loc:]
	pre = seq[loc-2:loc]
	if len(ss) > 4:
		subseq = ss[:4]
	else:
		subseq = ss
	presubseq = pre.lower()+subseq
	postsubseq = pre.lower() + 'N' + subseq[1] + 'T' + subseq[3:]

	ex = " "			
	for ch in presubseq:
		if hydropathy[ch.upper()] > 0:
			ex += "+"
		else:
			ex += "-"
	perfect = " ++-+-+"
	pattern_score = 0
	pattern_pos = 1
	while pattern_pos < len(perfect):
		if ex[pattern_pos] == perfect[pattern_pos]:
			pattern_score += 1
		pattern_pos += 1
		
	oldex = ex
	ex += " " +str(pattern_score)

	if oldex == " ++-+-+":
		ex += "  !!!"
	elif oldex == " -+-+-+":
		ex += "  ***"
	
	if presubseq[2] == "N" and (presubseq[4] == "S" or presubseq[4] == "T"):
		if 'P' not in presubseq.upper():
			ex += " already glycosylated"
	if loc in exclusions:
		ex += " (excluded)"
	elif loc in too_close:
		ex += " (too close)"

	return("; %5d %5d  %s --> %s %s\n" % (loc, val, presubseq, postsubseq, ex))
	
def dump_scores(scores, fh, seq, exclusions, too_close):
	fh.write(";\n;   loc score          mutation\n;\n")
			
	for score in scores:
		loc = score[0]
		val = score[1]
		str = format_score(loc, val, seq, exclusions, too_close)
		fh.write(str)
	
	fh.write(";\n")
	return
	
def dump_all_scores(dumpfile, scores, sorted_scores, seq, exclusions, too_close):
	with open(dumpfile, "w") as fh:
		dump_scores(scores, fh, seq, exclusions, too_close)
		dump_scores(sorted_scores, fh, seq, exclusions, too_close)
	return
	
def write_fasta(fn, scores, or_fasta_str, seq, exclusions, too_close):
	with open(fn, "w") as fh:
		fh.write(or_fasta_str)
		
		dump_scores(scores, fh, seq, exclusions, too_close)
		
		cnt = 0
		for ch in seq:
			cnt += 1
			fh.write(ch.upper())
			if cnt % 60 == 0:
				fh.write('\n')		
		if cnt % 60 != 0:
			fh.write('\n')
	return
			
def score_site(seq, n):
	
	if 'P' in seq[n-2:n+4]:
		return -999

	score = 0
	score += hydropathy[seq[n-2]]*10	# plus  - hydrophobic
	score += hydropathy[seq[n-1]]*10	# plus  - hydrophobic
	score -= hydropathy[seq[n]]*10		# minus - hydrophilic
	score += hydropathy[seq[n+1]]*10	# plus  - hydrophobic
	score -= hydropathy[seq[n+2]]*10	# minus - hydrophilic
	score += hydropathy[seq[n+3]]*10	# plus  - hydrophobic	
	return score
	
def glycosylate_protein(exclusions):
	scores = []
	top_sites = []
	too_close = []
	
# score sites
	
	fasta = read_fasta(infile)
	fasta_str = fasta[0]
	sequence = fasta[1]
	lower_limit = 2
	upper_limit = len(sequence) - 3
	n = lower_limit
	while n < upper_limit:
		site_score = score_site(sequence, n)
		scores.append([n, site_score])
		n += 1
		
	sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)

	good_sites = 0
	sites = 0
	for sc in sorted_scores:
		top = sc[0]
		for site in top_sites:
			if abs(top-site) < 3:
				too_close.append(top)
		top_sites.append(top)
		if top not in exclusions and top not in too_close:
			good_sites += 1
		if good_sites > num_sites: break
		sites += 1
		
	top_scores = sorted_scores[0:sites]
	
# generate mutant

	gly = ""

	n = 0
	while n < lower_limit:
		gly += sequence[n]
		n += 1

	n = lower_limit
	while n < upper_limit:
		if n in top_sites and n not in exclusions and n not in too_close:
			gly += 'n'	# asparagine
			n += 1
			gly += sequence[n]
			n += 1
			gly += 't'	# threonine
			n += 1
		else:
			gly += sequence[n]
			n += 1
			
	while n<len(sequence):
		gly += sequence[n]
		n += 1
			
	print(gly)
	print(len(gly))
	print(sequence)
	print(len(sequence))

	write_fasta(outfile, top_scores, fasta_str, gly, exclusions, too_close)
	too_close = []
	exclusions = []
	dump_all_scores(dumpfile, scores, sorted_scores, gly, exclusions, too_close)
	return
	
# main
	
glycosylate_protein(exclusions)	
