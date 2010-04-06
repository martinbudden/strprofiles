"""
strmarker
"""

from operator import itemgetter


SGM_PLUS_MARKERS = ['FGA','TH01','VWA','D2S1338','D3S1358','D8S1179','D16S539','D18S51','D19S433','D21S11']
#SGM_PLUS_MARKERS = ['TH01','D16S539']

def calculate_marker_rmp(alleles,theta):
	"""
	Calculate the random match probability for at a genetic marker, given the allele frequencies
	for that marker. The Balding and Nichols formulae are used.
	Keyword arguments:
	alleles -- a dict of allele frequencies in the form {'value':frequency,...}
	theta -- the population subdivision coefficient, used to correct for subdivided populations

	"""

	rmp = 0.0;
	denom = (1+theta)*(1+2*theta)
	for i in alleles:
		for j in alleles:
			if i==j:
				p = (2*theta+(1-theta)*alleles[i])*(3*theta+(1-theta)*alleles[i])/denom
				rmp += p*p
			else:
				p = (theta+(1-theta)*alleles[i])*(theta+(1-theta)*alleles[j])/denom
				rmp += 2*p*p
	return rmp


def pool_alleles(alleles,cutoff,count):
	"""
	Pool low frequency alleles together

	Keyword arguments:
	alleles -- a dict of allele frequencies in the form {'value':frequency,...}
	cutoff -- the minimum size of a frequency bin, items with a frequency lower than this will be pooled
	count -- the size of the sample from which the frequencies were derived

	"""

	items = alleles.items()
	items.sort(key = itemgetter(1))
	limit = float(cutoff)
	if count != 0:
		limit = float(cutoff)/float(count)
	otherCount = 0.0
	ret = {}
	for i in items:
		if i[1]<limit or otherCount<limit:
			otherCount += i[1]
		else:
			ret[i[0]] = i[1]
	if otherCount != 0:
		ret['other'] = otherCount
	return ret


def calc3(d,cols,name,cutoff,theta):
	rmp = 1.0
	col = {}
	for c in cols:
		#print c
		if c['name']==name and c['marker'] in SGM_PLUS_MARKERS:
			#print c['name'],c['marker'],c['count']
			count = c['count']
			# note, count doubled since two allele values per person in sample
			alleles = pool_alleles(c['alleles'],cutoff,2*count)
			p = calculate_marker_rmp(alleles,theta)
			rmp *= p
			#print c['name'],c['marker'],round(p,4)
			d[name][c['marker']] = p
	d[name]['count'] = count
	d[name]['combined'] = rmp
	d[name]['reciprocal'] = 1.0/rmp
	return d


def get_modal_profile(cols,name):
	profile = {}
	for c in cols:
		if c['name']==name and c['marker'] in SGM_PLUS_MARKERS:
			items = c['alleles'].items()
			items.sort(key = itemgetter(1))
			items.reverse()
			p = items[0][1]
			q = items[1][1]
			if p*p > 2*p*q:
				profile[c['marker']]  = ((items[0][0],p),(items[0][0],p))
			else:
				profile[c['marker']]  = ((items[0][0],p),(items[1][0],q))
	return profile


def calc_profile_match_probability(profile,theta):
	pmp = 1.0
	for marker in profile:
		m = profile[marker]
		p = m[0][1]
		q = m[1][1]
		#print m
		#mp *= 2*p*q
		denom = (1+theta)*(1+2*theta)
		if m[0][0] == m[1][0]:
			mp = (2*theta+(1-theta)*p)*(3*theta+(1-theta)*p)/denom
		else:
			mp = 2*(theta+(1-theta)*p)*(theta+(1-theta)*q)/denom
		pmp *= mp
	return pmp


