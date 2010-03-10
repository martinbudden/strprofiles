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

    def testRandomMatchProbability(self):
        """random match probability"""
	alleles = {'24': 0.1325, '25': 0.1125, '20.2': 0.0075, '27': 0.005, '20': 0.1625, '21': 0.1775, '22': 0.165, '23': 0.14, '19': 0.0625, '18': 0.015, '26': 0.015, '22.2': 0.005}
        self.assertTrue(abs(sgm.calculate_random_match_probability(alleles,0.0)-0.03557)<0.001, True)

    def testRandomMatchProbability2(self):
        """random match probability2"""
	pass

if __name__ == "__main__":
    unittest.main()