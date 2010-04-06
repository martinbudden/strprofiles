#!/usr/bin/env python
#coding=utf-8
#file: sgm.py

"""
sgm
"""

import array
import csv
import string
from operator import itemgetter
from jinja2 import Template
from collections import defaultdict
from optparse import OptionParser


SGM_PLUS_MARKERS = ['FGA','TH01','VWA','D2S1338','D3S1358','D8S1179','D16S539','D18S51','D19S433','D21S11']
#SGM_PLUS_MARKERS = ['TH01','D16S539']

def read_csv(filename,prefix,normalizer):
	"""
	Read genetic marker allele frequency data from a .csv file

	Keyword arguments:
	filename -- the file to be read
	normalizer -- used to normalize the data, use 100 if the data are percentages, 1 if they are
	fractions in the range [0.0,1.0]

	Format of file is:
	First row - each cell is a SGM Marker name, cell A1 blank
	Second row - each cell is sample size and sample name (eg "302 Cau"), cell B1 ignored
	First column - (starting at C1) allele value (eg "11"), cell A1 blank, cell B1 is column header "Allele"
	Main table - from C2 - allele frequency value for given Marker, Sample and allele

	,"CSF1PO","CSF1PO","CSF1PO","FGA","FGA","FGA",...
	"Allele","302 Cau","258 AA","140 His","302 Cau","258 AA","140 His",...
	7,,0.05253,0.02143,,,,...
	"""
	cols = []

	reader = csv.reader(open(filename, "rb"),'excel')

	# read in the first row which contains the
	row = reader.next()
	for i in range(1,len(row)-1):
		cols.append({'marker':row[i],'alleles':{}})
	row = reader.next()
	for i in range(1,len(row)-1):
		vals = string.split(row[i])
		if len(vals) == 1:
			cols[i-1]['name'] = prefix + vals[0]
		else:
			cols[i-1]['count'] = int(vals[0])
			cols[i-1]['name'] = prefix + vals[1]
	try:
		while 1:
			row = reader.next()
			for i in range(1,len(row)-1):
				if row[i] != '':
					cols[i-1]['alleles'][row[0]] = float(row[i])/normalizer
	except StopIteration:
		pass
	return cols


def textTable(data,caption,rowheaders,colheaders):
	"""
	Format data into a text table formatted using whitespace

	Keyword arguments:
	data -- a 2D dict of the data to be formatted
	caption -- the table's caption
	rowheaders -- the table's row headings
	colheaders -- the table's column headings

	"""
	t = Template('\n{{caption}}\n'
		'{{"%8s"|format("")}} '
		'{% for top in tops %}'
			'{{"%8s"|format(top)}} '
		'{% endfor %}\n'
		'{{"%8s"|format("Count")}} '
		'{% for top in tops %}'
			'{{"%8d"|format(data[top]["count"])|escape}} '
		'{% endfor %}\n'
		'{% for left in lefts %}'
			'{{"%8s"|format(left)}} '
			'{% for top in tops %}'
				'{{"   %1.3f"|format(data[top][left])|escape}} '
			'{% endfor %}\n'
		'{% endfor %}'
		'{{"%8s"|format("Combined")}} '
		'{% for top in tops %}'
			'{{"%1.2e"|format(data[top]["combined"])|escape}} '
		'{% endfor %}\n'
		'{{"%8s"|format("Recip")}} '
		'{% for top in tops %}'
			'{{"%1.2e"|format(data[top]["reciprocal"])|escape}} '
		'{% endfor %}')
	return t.render(data=data, caption=caption, lefts=rowheaders, tops=colheaders)


def textTable2(data,caption,rowheaders,colheaders):
	"""
	Format data into a text table formatted using whitespace

	Keyword arguments:
	data -- a 2D dict of the data to be formatted
	caption -- the table's caption
	rowheaders -- the table's row headings
	colheaders -- the table's column headings

	"""
	t = Template('\n{{caption}}\n'
		'{{"%8s"|format("")}} '
		'{% for top in tops %}'
			'{{"%8s"|format(top)}} '
		'{% endfor %}\n'
		'{% for left in lefts %}'
			'{{"%8s"|format(left)}} '
			'{% for top in tops %}'
				'{{"%1.2e"|format(data[top][left])|escape}} '
			'{% endfor %}\n'
		'{% endfor %}')
	return t.render(data=data, caption=caption, lefts=rowheaders, tops=colheaders)


def htmlTable(data,caption,rowheaders,colheaders):
	"""
	Format data into an HTML table

	Keyword arguments:
	data -- a 2D dict of the data to be formatted
	caption -- the table's caption
	rowheaders -- the table's row headings
	colheaders -- the table's column headings

	"""
	t = Template('<html>\n<table>\n'
		'<caption>{{caption}}</caption>\n'
		'<thead>\n'
		'<tr> '
		'<th></th> '
		'{% for top in tops %}'
			'<th>{{top|replace(" ", "<br />")}}</th> '
		'{% endfor %}'
		'</tr>\n'
		'<tr> '
		'<th>Count</th> '
		'{% for top in tops %}'
			'<th>{{data[top]["count"]|escape}}</th> '
		'{% endfor %}\n'
		'</tr>\n'
		'</thead>\n'
		'<tbody>\n'
		'{% for left in lefts %}'
			'<tr> '
			'<th>{{left}}</th> '
			'{% for top in tops %}'
				'<td>{{" %1.3f"|format(data[top][left])|escape}}</td> '
			'{% endfor %}'
			'</tr>\n'
		'{% endfor %}'
		'</tbody>\n'
		'<tr> '
		'<th>Combined</th> '
		'{% for top in tops %}'
			'<td>{{"%1.2e"|format(data[top]["combined"])|escape}}</td> '
		'{% endfor %}'
		'</tr>\n'
		'<tr> '
		'<th>Reciprocal</th> '
		'{% for top in tops %}'
			'<td>{{"%1.2e"|format(data[top]["reciprocal"])|escape}}</td> '
		'{% endfor %}'
		'</tr>\n'
		'<tfoot>\n'
		'</tfoot>\n'
		'</table>\n</html>\n')
	return t.render(data=data, caption=caption, lefts=rowheaders, tops=colheaders)


def htmlTable2(data,caption,rowheaders,colheaders):
	"""
	Format data into an HTML table

	Keyword arguments:
	data -- a 2D dict of the data to be formatted
	caption -- the table's caption
	rowheaders -- the table's row headings
	colheaders -- the table's column headings

	"""
	t = Template('<html>\n<table>\n'
		'<caption>{{caption}}</caption>\n'
		'<thead>\n'
		'<tr> '
		'<th></th> '
		'{% for top in tops %}'
			'<th>{{top|replace(" ", "<br />")}}</th> '
		'{% endfor %}'
		'</tr>\n'
		'</thead>\n'
		'<tbody>\n'
		'{% for left in lefts %}'
			'<tr> '
			'<th>{{left}}</th> '
			'{% for top in tops %}'
				'<td>{{"%1.2e"|format(data[top][left])|escape}}</td> '
			'{% endfor %}'
			'</tr>\n'
		'{% endfor %}'
		'</tbody>\n'
		'</table>\n</html>\n')
	return t.render(data=data, caption=caption, lefts=rowheaders, tops=colheaders)


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


def calc4(cols,samples,format,caption,cutoff,theta):
	d = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		calc3(d,cols,sample,cutoff,theta)
	#"302 Cau","258 AA","140 His",
	#tops = ['AB AA','AB Cau','AB Combined','JSP AA','JSP Cau','JSP His','JSP Combined']
	#columnHeaders = ['AB AA','AB Cau','JSF AA','JSF Cau','JSF His']
	if format == "html":
		return htmlTable(d,caption,SGM_PLUS_MARKERS,columnHeaders)
	else:
		return textTable(d,caption,SGM_PLUS_MARKERS,columnHeaders)


def calc5(cols,samples,format,caption):
	d = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		profile = get_modal_profile(cols,sample)
		d[sample]['0.0'] = 1.0/calc_profile_match_probability(profile,0.0)
		d[sample]['0.01'] = 1.0/calc_profile_match_probability(profile,0.01)
		d[sample]['0.03'] = 1.0/calc_profile_match_probability(profile,0.03)
	if format == "html":
		return htmlTable2(d,caption,['0.0','0.01','0.03'],columnHeaders)
	else:
		return textTable2(d,caption,['0.0','0.01','0.03'],columnHeaders)


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


def main():
	parser = OptionParser()
	parser.add_option("-t", action="store_true", dest="text_format", default=False, help="use text format tables")
 	parser.add_option("-v", action="store_true", dest="verbose", default=False, help="print status messages to stdout")
	(options,args) = parser.parse_args()
	if options.text_format:
		format = "text"
	else:
		format = "html"

	# read in the NIST/JSF allele frequency data
	colsJFS = read_csv("../data/JFS2003IDresults.csv","JSF ",1)

	# read in the NIST/JSF allele frequency data
	colsAB = read_csv("../data/ABresults.csv","AB ",100)
	cols = colsJFS + colsAB

	samples = ['JSF AA','JSF Cau','JSF His','AB AA','AB Cau']

	print calc4(cols,samples,format,"Raw Probability of Identity values",0,0.0)
	print calc4(cols,samples,format,"Rare alleles pooled",5,0.0)
	print calc4(cols,samples,format,"Theta = 0.01",5,0.01)
	print calc4(cols,samples,format,"Theta = 0.03",5,0.03)
	print calc5(cols,samples,format,"Modal Man")

if __name__ == "__main__":
	main()