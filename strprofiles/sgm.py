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
	data = []
	reader = csv.reader(open(filename, "rb"),'excel')

	# read in the first row, which contains the SMG Marker names
	row = reader.next()
	for i in range(1,len(row)-1):
		data.append({'marker':row[i],'alleles':{}})

	# read in the second row, which contains the sample names
	row = reader.next()
	for i in range(1,len(row)-1):
		vals = string.split(row[i])
		if len(vals) == 1:
			data[i-1]['name'] = prefix + vals[0]
		else:
			data[i-1]['count'] = int(vals[0])
			data[i-1]['name'] = prefix + vals[1]

	# read in the allele frequencies
	try:
		while 1:
			row = reader.next()
			for i in range(1,len(row)-1):
				if row[i] != '':
					data[i-1]['alleles'][row[0]] = float(row[i])/normalizer
	except StopIteration:
		pass
	return data


def textRmpTable(data,caption,rowheaders,colheaders):
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
		'{% for colheader in colheaders %}'
			'{{"%8s"|format(colheader)}} '
		'{% endfor %}\n'
		'{{"%8s"|format("Count")}} '
		'{% for colheader in colheaders %}'
			'{{"%8d"|format(data[colheader]["count"])|escape}} '
		'{% endfor %}\n'
		'{% for rowheader in rowheaders %}'
			'{{"%8s"|format(rowheader)}} '
			'{% for colheader in colheaders %}'
				'{{"   %1.3f"|format(data[colheader][rowheader])|escape}} '
			'{% endfor %}\n'
		'{% endfor %}'
		'{{"%8s"|format("Combined")}} '
		'{% for colheader in colheaders %}'
			'{{"%1.2e"|format(data[colheader]["combined"])|escape}} '
		'{% endfor %}\n'
		'{{"%8s"|format("Recip")}} '
		'{% for colheader in colheaders %}'
			'{{"%1.2e"|format(data[colheader]["reciprocal"])|escape}} '
		'{% endfor %}')
	return t.render(data=data, caption=caption, rowheaders=rowheaders, colheaders=colheaders)


def textPmpTable(data,caption,rowheaders,colheaders):
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
		'{% for colheader in colheaders %}'
			'{{"%8s"|format(colheader)}} '
		'{% endfor %}\n'
		'{% for rowheader in rowheaders %}'
			'{{"%8s"|format(rowheader)}} '
			'{% for colheader in colheaders %}'
				'{{"%1.2e"|format(data[colheader][rowheader])|escape}} '
			'{% endfor %}\n'
		'{% endfor %}')
	return t.render(data=data, caption=caption, rowheaders=rowheaders, colheaders=colheaders)


def htmlRmpTable(data,caption,rowheaders,colheaders):
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
		'{% for colheader in colheaders %}'
			'<th>{{colheader|replace(" ", "<br />")}}</th> '
		'{% endfor %}'
		'</tr>\n'
		'<tr> '
		'<th>Count</th> '
		'{% for colheader in colheaders %}'
			'<th>{{data[colheader]["count"]|escape}}</th> '
		'{% endfor %}\n'
		'</tr>\n'
		'</thead>\n'
		'<tbody>\n'
		'{% for rowheader in rowheaders %}'
			'<tr> '
			'<th>{{rowheader}}</th> '
			'{% for colheader in colheaders %}'
				'<td>{{" %1.3f"|format(data[colheader][rowheader])|escape}}</td> '
			'{% endfor %}'
			'</tr>\n'
		'{% endfor %}'
		'</tbody>\n'
		'<tr> '
		'<th>Combined</th> '
		'{% for colheader in colheaders %}'
			'<td>{{"%1.2e"|format(data[colheader]["combined"])|escape}}</td> '
		'{% endfor %}'
		'</tr>\n'
		'<tr> '
		'<th>Reciprocal</th> '
		'{% for colheader in colheaders %}'
			'<td>{{"%1.2e"|format(data[colheader]["reciprocal"])|escape}}</td> '
		'{% endfor %}'
		'</tr>\n'
		'<tfoot>\n'
		'</tfoot>\n'
		'</table>\n</html>\n')
	return t.render(data=data, caption=caption, rowheaders=rowheaders, colheaders=colheaders)


def htmlPmpTable(data,caption,rowheaders,colheaders):
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
		'{% for colheader in colheaders %}'
			'<th>{{colheader|replace(" ", "<br />")}}</th> '
		'{% endfor %}'
		'</tr>\n'
		'</thead>\n'
		'<tbody>\n'
		'{% for rowheader in rowheaders %}'
			'<tr> '
			'<th>{{rowheader}}</th> '
			'{% for colheader in colheaders %}'
				'<td>{{"%1.2e"|format(data[colheader][rowheader])|escape}}</td> '
			'{% endfor %}'
			'</tr>\n'
		'{% endfor %}'
		'</tbody>\n'
		'</table>\n</html>\n')
	return t.render(data=data, caption=caption, rowheaders=rowheaders, colheaders=colheaders)


def calc_rmps(data,samples,format,caption,cutoff,theta):
	"""
	Calculate the random match probabilities and format them into a table
	"""
	table = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		table[sample] = strmarker.calc_rmps(data,sample,cutoff,theta)
	if format == "html":
		return htmlRmpTable(table,caption,strmarker.SGM_PLUS_MARKERS,columnHeaders)
	else:
		return textRmpTable(table,caption,strmarker.SGM_PLUS_MARKERS,columnHeaders)


def calc_pmps(data,samples,format,caption):
	"""
	Calculate the profile match probabilities for the modal profile and format them into a table
	"""
	table = defaultdict(dict)
	columnHeaders = []
	for sample in samples:
		columnHeaders.append(sample)
		profile = strmarker.get_modal_profile(data,sample)
		table[sample]['0.0'] = 1.0/strmarker.calc_profile_match_probability(profile,0.0)
		table[sample]['0.01'] = 1.0/strmarker.calc_profile_match_probability(profile,0.01)
		table[sample]['0.03'] = 1.0/strmarker.calc_profile_match_probability(profile,0.03)
	if format == "html":
		return htmlPmpTable(table,caption,['0.0','0.01','0.03'],columnHeaders)
	else:
		return textPmpTable(table,caption,['0.0','0.01','0.03'],columnHeaders)


def main():
	"""
	Read in the allele frequency data
	Print out tables of random match probabilities.
	Print out the table of profile match probabilities for the modal profile.

	"""

	parser = OptionParser()
	parser.add_option("-t", action="store_true", dest="text_format", default=False, help="use text format tables")
 	parser.add_option("-v", action="store_true", dest="verbose", default=False, help="print status messages to stdout")
	(options,args) = parser.parse_args()
	if options.text_format:
		format = "text"
	else:
		format = "html"

	# read in the NIST/JSF allele frequency data
	dataJFS = read_csv("../data/JFS2003IDresults.csv","JSF ",1)

	# read in the NIST/JSF allele frequency data
	dataAB = read_csv("../data/ABresults.csv","AB ",100)
	data = dataJFS + dataAB

	samples = ['JSF AA','JSF Cau','JSF His','AB AA','AB Cau']

	print calc_rmps(data,samples,format,"Raw Probability of Identity values",0,0.0)
	print calc_rmps(data,samples,format,"Rare alleles pooled",5,0.0)
	print calc_rmps(data,samples,format,"Theta = 0.01",5,0.01)
	print calc_rmps(data,samples,format,"Theta = 0.03",5,0.03)
	print calc_pmps(data,samples,format,"Modal Man")


if __name__ == "__main__":
	main()