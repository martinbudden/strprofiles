"""
Example test module.
"""

import unittest

from strprofiles import sgm

class exampleTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testMarkerRandomMatchProbability(self):
        """marker random match probability"""
	# AB Cau TH01 frequencies (n=200)
	alleles = self.fromPercentages({'6':25.25, '7':16.00, '8':9.25, '9':13.75, '9.3':35.00, 10:0.75})
	expected = 0.094
	result = sgm.calculate_marker_rmp(alleles,0.0)
	print result
        self.assertTrue(round(result,3)==expected, True)

	# AB Cau FGA frequencies (n=200)
	alleles = {'18': 0.015, '19': 0.0625, '20': 0.1625, '20.2': 0.0075, '21': 0.1775, '22': 0.165, '22.2': 0.005, '23': 0.14, '24': 0.1325, '25': 0.1125, '26': 0.015, '27': 0.005}
	expected = 0.036
	result = sgm.calculate_marker_rmp(alleles,0.0)
        self.assertTrue(round(result,3)==expected, True)

	# AB Cau D16S539 frequencies (n=200)
	alleles = self.fromPercentages({'8':1.5,'9':10.75,'10':4.25,'11':29.75,'12':34.25,'13':17.5,'14':1.75,'15':0.25})
	expected = 0.103
	result = sgm.calculate_marker_rmp(alleles,0.0)
        self.assertTrue(round(result,3)==expected, True)

    def fromPercentages(self,alleles):
        for i in alleles:
            alleles[i] /= 100.0
        return alleles

    def testPoolAlleles(self):
        """pool alleles"""
	# constructed counts (n=100)
	alleles = {'a':1,'b':2,'c':3,'d':4,'e':5,'f':15,'g':30,'h':40}
	self.assertTrue(sgm.pool_alleles(alleles,0,0)==alleles)
	expected = {'other':10,'e':5,'f':15,'g':30,'h':40}
	self.assertTrue(sgm.pool_alleles(alleles,5,0)==expected)

	# constructed frequencies (n=100)
	alleles = {'a':0.01,'b':0.02,'c':0.03,'d':0.04,'e':0.05,'f':0.15,'g':0.30,'h':0.40}
	result = sgm.pool_alleles(alleles,0,100)
	for i in result:
		self.assertTrue(abs(result[i]-alleles[i])<0.0001)
	expected = {'other':0.10,'e':0.05,'f':0.15,'g':0.30,'h':0.40}
	result = sgm.pool_alleles(alleles,5,100)
	for i in result:
		self.assertTrue(abs(result[i]-expected[i])<0.0001)

	# AB Cau D16S539 counts (n=200)
	alleles = {'8':6,'9':43,'10':17,'11':119,'12':137,'13':70,'14':7,'15':1}
	self.assertTrue(sgm.pool_alleles(alleles,0,0)==alleles)
	expected = {'9':43,'10':17,'11':119,'12':137,'13':70,'14':7,'other':7}
	self.assertTrue(sgm.pool_alleles(alleles,5,0)==expected)

	# AB Cau D16S539 frequencies (n=200)
	alleles = self.fromPercentages({'8':1.5,'9':10.75,'10':4.25,'11':29.75,'12':34.25,'13':17.5,'14':1.75,'15':0.25})
	result = sgm.pool_alleles(alleles,0,200)
	for i in result:
		self.assertTrue(abs(result[i]-alleles[i])<0.0001)
	expected = self.fromPercentages({'9':10.75,'10':4.25,'11':29.75,'12':34.25,'13':17.5,'14':1.75,'other':1.75})
	result = sgm.pool_alleles(alleles,5,200*2)
	for i in result:
		self.assertTrue(abs(result[i]-expected[i])<0.0001)



if __name__ == "__main__":
    unittest.main()