"""
strmarker

Utilities for STR (short tandem repeat) genetic markers and allele frequency data.
"""

from operator import itemgetter


SGM_PLUS_MARKERS = ['FGA','TH01','VWA','D2S1338','D3S1358','D8S1179','D16S539','D18S51','D19S433','D21S11']

def calc_marker_rmp(alleles,theta):
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


def calc_rmps(data,name,cutoff,theta):
	"""
	Calculate the random match probabilities for each genetic marker for the named sample set

	Keyword arguments:
	data -- the allele frequency data
	name -- the name of the sample set to be used
	cutoff -- the minimum size of a frequency bin, items with a frequency lower than this will be pooled
	theta -- the population subdivision coefficient, used to correct for subdivided populations

	"""
	ret = {}
	rmp = 1.0
	for d in data:
		if d['name']==name and d['marker'] in SGM_PLUS_MARKERS:
			count = d['count']
			# note, count doubled since two allele values per person in sample
			alleles = pool_alleles(d['alleles'],cutoff,2*count)
			p = calc_marker_rmp(alleles,theta)
			ret[d['marker']] = p
			rmp *= p
	ret['count'] = count
	ret['combined'] = rmp
	ret['reciprocal'] = 1.0/rmp
	return ret


def get_modal_profile(data,name):
	"""
	Find the modal profile for in the named sample set.

	"""
	profile = {}
	for d in data:
		if d['name']==name and d['marker'] in SGM_PLUS_MARKERS:
			items = d['alleles'].items()
			items.sort(key = itemgetter(1))
			items.reverse()
			p = items[0][1]
			q = items[1][1]
			if p*p > 2*p*q:
				profile[d['marker']]  = ((items[0][0],p),(items[0][0],p))
			else:
				profile[d['marker']]  = ((items[0][0],p),(items[1][0],q))
	return profile


def calc_profile_match_probability(profile,theta):
	"""
	Calculate the probability that a random individual matches the given profile

	Keyword arguments:
	profile -- a 2D dict of the data to be formatted
	theta -- the population subdivision coefficient, used to correct for subdivided populations

	"""
	pmp = 1.0
	for marker in profile:
		m = profile[marker]
		p = m[0][1]
		q = m[1][1]
		denom = (1+theta)*(1+2*theta)
		if m[0][0] == m[1][0]:
			mp = (2*theta+(1-theta)*p)*(3*theta+(1-theta)*p)/denom
		else:
			mp = 2*(theta+(1-theta)*p)*(theta+(1-theta)*q)/denom
		pmp *= mp
		#print marker,round(mp,4)
	return pmp


