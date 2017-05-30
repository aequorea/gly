# gly1.py -- protein glycosylator
# 2016-10-05 -- first version (John Saeger)

# This version attempts to implement what is known about glycosylation
# sites as reported in this study: 

# The Structural Assessment of Glycosylation Sites Database - SAGS -
# An Overall View on N-Glycosylation by Marius D. Surleac, 
# Laurentiu N. Spiridon, Robi Tacutu, Adina L Milac, 
# Stefana M. Petrescu and Andrei-J Petrescu 
# (http://dx.doi.org/10.5772/51690)

# The scoring algorithm assigns points based on how many standard 
# deviations away from normal certain types of residues occur 
# in the neighborhood of a glycosylation site.  
# We give a strong preference to neutral polar residues at 
# position n and n+2. We exclude prolines from n-2 to n+3.


infile = "eyfp.fas"
outfile = "eyfp.gly1"		# first version output
num_sites = 20				# desired number of glycosylation sites
exclusions = [247]			# mutation sites to disallow

hydrophobic = set('VLIM')
aromatic = set('FYW')		# Y - aromaticity overrides polarity?
acidic = set('ED')
basic = set('RHK')
neutral = set('NQST')

scores = []
top_sites = []
top_scores = []
gly = ""

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
	
def format_score(loc, val, seq):
	ss = sequence[loc:]
	if loc > 1:
		pre = sequence[loc-2:loc]
	elif loc > 0:
		pre = sequence[loc-1:loc]
	else:
		pre = ""
	if len(ss) > 8:
		subseq = ss[:8]
	else:
		subseq = ss
	presubseq = pre.lower()+subseq
	postsubseq = pre.lower() + 'N' + subseq[1] + 'T' + subseq[3:]
	if loc in exclusions:
		ex = "(excluded)"
	else:
		ex = ""
	return("; %5d %5d  %s --> %s %s\n" % (loc, val, presubseq, postsubseq, ex))
	
def write_fasta(fn, scores, or_fasta_str, seq):
	with open(fn, "w") as fh:
		fh.write(or_fasta_str)
		fh.write(";\n;   loc score          mutation\n;\n")
				
		for score in scores:
			loc = score[0]
			val = score[1]
			str = format_score(loc, val, seq)
			fh.write(str)
		
		fh.write(";\n")

		cnt = 0
		for ch in seq:
			cnt += 1
			fh.write(ch.upper())
			if cnt % 60 == 0:
				fh.write('\n')		
		if cnt % 60 != 0:
			fh.write('\n')
			
def score_site(seq, n):		# could be a little more elegant
	score = 0
	if len(sequence) > n+3:
		if seq[n+3] in hydrophobic:
			score += 20
		elif seq[n+3] in aromatic:
			score += 15
		elif seq[n+3] in acidic|basic:
			score -= 5
		elif seq[n+3] == 'P':
			score = -999
		
	if len(sequence) > n+2:
		if seq[n+2] in neutral:
			score += 40
		elif seq[n+2] == 'P':
			score = -999
			
	if len(sequence) > n+1:
		if seq[n+1] in hydrophobic:
			score += 25
		elif seq[n+1] in aromatic:
			score += 5
		elif seq[n+1] in acidic|basic:
			score -= 15
		elif seq[n+1] == 'P':
			score = -999
		
	if seq[n] in neutral:
		score += 40
	elif seq[n] == 'P':
		score = -999
		
	if n>1:
		if seq[n-1] in hydrophobic:
			score += 8
		elif seq[n-1] in aromatic:
			score += 20
		elif seq[n-1] in acidic|basic:
			score -= 22
		elif seq[n-1] == 'P':
			score = -999
			
	if n>2:
		if seq[n-2] in hydrophobic:
			score += 5
		elif seq[n-2] in aromatic:
			score += 28
		elif seq[n-2] in acidic|basic:
			score -= 23
		elif seq[n-2] == 'P':
			score = -999	
	
	return score
	

ungly = read_fasta(infile)
fasta_str = ungly[0]
sequence = ungly[1]

n = 0
for ch in sequence:
	site_score = score_site(sequence, n)
	scores.append([n, site_score])
	n += 1
	
sorted_scores = sorted(scores, key=lambda s: s[1], reverse=True)

good_sites = 0
sites = 0
for sc in sorted_scores:
	top = sc[0]
	top_sites.append(top)
	if top not in exclusions:
		good_sites += 1
	if good_sites > num_sites: break
	sites += 1
	
top_scores = sorted_scores[0:sites]
	
pos = 0
while pos < len(sequence):
	if pos in top_sites and pos not in exclusions:
		gly += 'n'	# asparagine
		pos += 1
		gly += sequence[pos]
		pos += 1
		gly += 't'	# threonine
		pos += 1
	else:
		gly += sequence[pos]
		pos += 1
		
print(gly)

write_fasta(outfile, top_scores, fasta_str, gly)
