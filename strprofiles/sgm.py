#!/usr/bin/env python
#coding=utf-8
#file: sgm.py

"""
sgm
"""

import strmarker
import csv
import string
from jinja2 import Template
from collections import defaultdict
from optparse import OptionParser


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


def calc4(cols,samples,format,caption,cutoff,theta):
	d = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		d[sample] = strmarker.calc_rmps(cols,sample,cutoff,theta)
	#"302 Cau","258 AA","140 His",
	#tops = ['AB AA','AB Cau','AB Combined','JSP AA','JSP Cau','JSP His','JSP Combined']
	#columnHeaders = ['AB AA','AB Cau','JSF AA','JSF Cau','JSF His']
	if format == "html":
		return htmlTable(d,caption,strmarker.SGM_PLUS_MARKERS,columnHeaders)
	else:
		return textTable(d,caption,strmarker.SGM_PLUS_MARKERS,columnHeaders)


def calc5(cols,samples,format,caption):
	d = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		profile = strmarker.get_modal_profile(cols,sample)
		d[sample]['0.0'] = 1.0/strmarker.calc_profile_match_probability(profile,0.0)
		d[sample]['0.01'] = 1.0/strmarker.calc_profile_match_probability(profile,0.01)
		d[sample]['0.03'] = 1.0/strmarker.calc_profile_match_probability(profile,0.03)
	if format == "html":
		return htmlTable2(d,caption,['0.0','0.01','0.03'],columnHeaders)
	else:
		return textTable2(d,caption,['0.0','0.01','0.03'],columnHeaders)


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