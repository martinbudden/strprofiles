"""
sgm
"""

import array
import csv
from operator import itemgetter

SGM_PLUS_MARKERS = ['FGA','TH01','VWA','D2S1338','D3S1358','D8S1179','D16S539','D18S51','D19S433','D21S11']

def calculate_random_match_probability(alleles,theta):
	"""
	Calculate the random match probability given the allele frequencies
	"""
	rmp = 0.0;
	denom = (1+theta)*(1+2*theta)
	for i in alleles:
		a = alleles[i]
		for j in alleles:
			if i==j:
				p = (2*theta+(1-theta)*alleles[i])*(3*theta+(1-theta)*alleles[i])/denom
				rmp += p*p
			else:
				p = (theta+(1-theta)*alleles[i])*(theta+(1-theta)*alleles[j])/denom
				rmp += 2*p*p
	return rmp

def read_csv(filename,divisor):
	reader = csv.reader(open(filename, "rb"),'excel')
	row = reader.next()
	cols = []
	for i in range(1,len(row)-1):
		cols.append({'marker':row[i],'alleles':{}})
	row = reader.next()
	for i in range(1,len(row)-1):
		cols[i-1]['name']=row[i]
	try:
		while 1:
			row = reader.next()
			for i in range(1,len(row)-1):
				if row[i] != '':
					cols[i-1]['alleles'][row[0]] = float(row[i])/divisor
	except StopIteration:
		pass
	return cols

def groupAlleles(alleles,cutoff):
	items = alleles.items()
	items.sort(key = itemgetter(1))
	bin = 0.0
	ret = {}
	for i in items:
		if i[1]<cutoff or bin<cutoff:
			bin += i[1]
		else:
			ret[i[0]] = i[1]
	if bin!=0.0:
		ret['bin'] = bin
	#for i in ret:
	#	print i,':',round(ret[i],4),' ',
	#print
	return ret

def calc3(cols,name,cutoff,theta):
	rmp = 1.0
	for c in cols:
		if c['name']==name and c['marker'] in SGM_PLUS_MARKERS:
			alleles = groupAlleles(c['alleles'],cutoff)
			p = calculate_random_match_probability(alleles,theta)
			rmp *= p
			print c['name'],c['marker'],round(p,4)
	print rmp
	print "%.4e" % (1.0/rmp)


def main():
	#cols = read_csv("../data/JFS2003IDresults.csv",1)
	cols = read_csv("../data/ABresults.csv",100)
	#calc3(cols,'AA')
	#theta = 0.01
	#calc3(cols,'Cau',5.0/195.0,theta)
	calc3(cols,'Cau',0.0,0.0)

main()
