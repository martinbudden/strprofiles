"""
Example test module.
"""

import unittest

from strprofiles import strmarker

class exampleTestCase(unittest.TestCase):

    def setUp(self):
	# AB Cau FGA frequencies (n=200)
	self.AB_Cau_FGA_alleles = {'18': 0.015, '19': 0.0625, '20': 0.1625, '20.2': 0.0075, '21': 0.1775, '22': 0.165, '22.2': 0.005, '23': 0.14, '24': 0.1325, '25': 0.1125, '26': 0.015, '27': 0.005}
	# AB Cau TH01 frequencies (n=200)
	self.AB_Cau_TH01_alleles = self.fromPercentages({'6':25.25, '7':16.00, '8':9.25, '9':13.75, '9.3':35.00, 10:0.75})
	# AB Cau D16S539 frequencies (n=200)
	self.AB_Cau_D16S539_alleles = self.fromPercentages({'8':1.5,'9':10.75,'10':4.25,'11':29.75,'12':34.25,'13':17.5,'14':1.75,'15':0.25})


    def tearDown(self):
        pass

    def testMarkerRandomMatchProbability(self):
        """marker random match probability"""

	# AB Cau FGA frequencies (n=200)
	alleles = self.AB_Cau_FGA_alleles
	expected = 0.036
	result = strmarker.calc_marker_rmp(alleles,0.0)
        self.assertAlmostEqual(result,expected,3)

	# AB Cau TH01 frequencies (n=200)
	alleles = self.AB_Cau_TH01_alleles
	expected = 0.094
	result = strmarker.calc_marker_rmp(alleles,0.0)
        self.assertAlmostEqual(result,expected,3)

	# AB Cau D16S539 frequencies (n=200)
	alleles = self.AB_Cau_D16S539_alleles
	expected = 0.103
	result = strmarker.calc_marker_rmp(alleles,0.0)
        self.assertAlmostEqual(result,expected,3)


    def fromPercentages(self,alleles):
        for i in alleles:
            alleles[i] /= 100.0
        return alleles


    def testGetModalProfile(self):
        """get modal profile"""
	cols = [{'name':'AB','marker':'FGA','alleles':self.AB_Cau_FGA_alleles},{'name':'AB','marker':'TH01','alleles':self.AB_Cau_TH01_alleles}, \
		{'name':'AB','marker':'DS16S539','alleles':self.AB_Cau_D16S539_alleles}]
	expected = {'FGA':(('21',0.1775),('22',0.165)), 'TH01':(('9.3',0.35),('6',0.2525)),' DS16S539':(('12',0.3425),('11',0.2975))}
	profile = expected
	result = strmarker.get_modal_profile(cols,'AB')
	for i in result:
		self.assertEqual(result[i][0][0],expected[i][0][0])
		self.assertEqual(result[i][1][0],expected[i][1][0])
		self.assertAlmostEqual(result[i][0][1],expected[i][0][1])
		self.assertAlmostEqual(result[i][1][1],expected[i][1][1])
	expected = 0.01035313125
	result = strmarker.calc_profile_match_probability(result,0.0)
	self.assertEqual(result,expected)

    def testPoolAlleles(self):
        """pool alleles"""
	# constructed counts (n=100)
	alleles = {'a':1,'b':2,'c':3,'d':4,'e':5,'f':15,'g':30,'h':40}
	self.assertEqual(strmarker.pool_alleles(alleles,0,0),alleles)
	expected = {'other':10,'e':5,'f':15,'g':30,'h':40}
	self.assertEqual(strmarker.pool_alleles(alleles,5,0),expected)

	# constructed frequencies (n=100)
	alleles = {'a':0.01,'b':0.02,'c':0.03,'d':0.04,'e':0.05,'f':0.15,'g':0.30,'h':0.40}
	result = strmarker.pool_alleles(alleles,0,100)
	for i in result:
		self.assertAlmostEqual(result[i],alleles[i])
	expected = {'other':0.10,'e':0.05,'f':0.15,'g':0.30,'h':0.40}
	result = strmarker.pool_alleles(alleles,5,100)
	for i in result:
		self.assertAlmostEqual(result[i],expected[i])

	# AB Cau D16S539 counts (n=200)
	alleles = {'8':6,'9':43,'10':17,'11':119,'12':137,'13':70,'14':7,'15':1}
	self.assertEqual(strmarker.pool_alleles(alleles,0,0),alleles)
	expected = {'9':43,'10':17,'11':119,'12':137,'13':70,'14':7,'other':7}
	self.assertEqual(strmarker.pool_alleles(alleles,5,0),expected)

	# AB Cau D16S539 frequencies (n=200)
	alleles = self.AB_Cau_D16S539_alleles
	result = strmarker.pool_alleles(alleles,0,200)
	for i in result:
		self.assertAlmostEqual(result[i],alleles[i])
	expected = self.fromPercentages({'9':10.75,'10':4.25,'11':29.75,'12':34.25,'13':17.5,'14':1.75,'other':1.75})
	result = strmarker.pool_alleles(alleles,5,200*2)
	for i in result:
		self.assertAlmostEqual(result[i],expected[i])



if __name__ == "__main__":
		unittest.main()