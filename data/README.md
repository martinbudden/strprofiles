STR Genetic Marker Allele Frequency Data
========================================

There are two sets of data in this directory.


Data Format
===========

Format of file is:
First row - each cell is a SGM Marker name, cell A1 blank
Second row - each cell is sample size and sample name (eg "302 Cau"), cell B1 ignored
First column - (starting at cell C1) allele value (eg "11"), cell A1 blank, cell B1 is column header "Allele"
Main table - from cell C2 onwards - allele frequency value for given Marker, Sample and allele

        , "CSF1PO",  "CSF1PO", "CSF1PO",  "FGA",     "FGA",    "FGA",...
"Allele", "302 Cau", "258 AA", "140 His", "302 Cau", "258 AA", "140 His",...
7,                 ,  0.05253,   0.02143,...
...

NIST Data
=========

File: JFS2003IDresults.csv

This dataset is from the [Population Studies Conducted by the [National Institute of Standards and Technology (NIST) Forensics/Human Identity Project Team](http://www.cstl.nist.gov/biotech/strbase/NISTpop.htm). The team derived their data from anonymous male samples from the Interstate Blood Bank (Memphis, TN). These samples came from 258 African American,  302 U.S. Caucasian, and 140 Hispanic males (self-identified), see [Allele Frequencies for 15 Autosomal STR Loci on U.S. Caucasian, African American, and Hispanic Population](http://www.cstl.nist.gov/biotech/strbase/pub_pres/Butler2003a.pdf) The raw allele frequency data are available as an [Excel file](http://www.cstl.nist.gov/biotech/strbase/NISTpopdata/JFS2003IDresults.xls).

JFS2003IDresults.csv is the file JFS2003IDresults.xls saved in CSV format.

I refer to this data with the abbreviations "JFS AA", "JFS Cau" and "JSF His".


Applied Biosystems Data
=======================

File: ABresults.csv

This dataset is that used by Applied Biosystems in their [AmpFlSTR® SGM Plus® PCR Amplification Kit User’s Manual](http://www3.appliedbiosystems.com/cms/groups/applied_markets_support/documents/generaldocuments/cms_041049.pdf). These data were derived from samples provided by the Laboratory Corporation of America. The samples came from 195 African American and 200 U.S. Caucasian individuals throughout the United States with no geographical preference.

I have manually transcribed this data and saved in the file ABresults.csv, in the same format as the NIST Data.

I refer to this data with the abbreviations "AB AA" and "AB Cau".

