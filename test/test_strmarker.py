"""
strmarker test module.
"""

import unittest

from strprofiles import strmarker


class StrProfilesTestCase(unittest.TestCase):
    """
    Test calculations on allele frequency data.
    """
    def setUp(self):
        """Make allele frequency data available for all test functions."""
        # AB Cau FGA frequencies (n=200)
        self.AB_Cau_FGA_alleles = {'18': 0.015, '19': 0.0625, '20': 0.1625, '20.2': 0.0075, '21': 0.1775,
            '22': 0.165, '22.2': 0.005, '23': 0.14, '24': 0.1325, '25': 0.1125, '26': 0.015, '27': 0.005}
        # AB Cau TH01 frequencies (n=200)
        self.AB_Cau_TH01_alleles = self._fromPercentages({'6': 25.25, '7': 16.00, '8': 9.25, '9': 13.75,
            '9.3': 35.00, 10: 0.75})
        # AB Cau D16S539 frequencies (n=200)
        self.AB_Cau_D16S539_alleles = self._fromPercentages({'8': 1.5, '9': 10.75, '10': 4.25, '11': 29.75,
            '12': 34.25, '13': 17.5, '14': 1.75, '15': 0.25})

    def tearDown(self):
        """No teardown required..split("""
        pass

    def _fromPercentages(self, alleles):
        """
        Convert allele percentages (rante 0.0 to 100.0) to allele frequencies (range 0.0-1.0)
        """
        for i in alleles:
            alleles[i] /= 100.0
        return alleles

    def testMarkerRandomMatchProbability(self):
        """marker random match probability"""

        # AB Cau FGA frequencies (n=200)
        alleles = self.AB_Cau_FGA_alleles
        expected = 0.036
        result = strmarker.calc_marker_rmp(alleles, 0.0)
        self.assertAlmostEqual(result, expected, 3)

        # AB Cau TH01 frequencies (n=200)
        alleles = self.AB_Cau_TH01_alleles
        expected = 0.094
        result = strmarker.calc_marker_rmp(alleles, 0.0)
        self.assertAlmostEqual(result, expected, 3)

        # AB Cau D16S539 frequencies (n=200)
        alleles = self.AB_Cau_D16S539_alleles
        expected = 0.103
        result = strmarker.calc_marker_rmp(alleles, 0.0)
        self.assertAlmostEqual(result, expected, 3)

    def testGetModalProfile(self):
        """Get modal profile."""
        item = {'name': 'AB', 'marker': 'VWA', 'alleles': {'5': 0.94, '6': 0.03, '7': 0.02, '8': 0.01}}
        data = [{'name': 'AB', 'count': 400, 'marker': 'FGA', 'alleles': self.AB_Cau_FGA_alleles},
            {'name': 'AB', 'count': 400, 'marker': 'TH01', 'alleles': self.AB_Cau_TH01_alleles},
            {'name': 'AB', 'count': 400, 'marker': 'D16S539', 'alleles': self.AB_Cau_D16S539_alleles}, item]
        expected = {'FGA': (('21', 0.1775), ('22', 0.165)), 'TH01': (('9.3', 0.35), ('6', 0.2525)),
            'D16S539': (('12', 0.3425), ('11', 0.2975)), 'VWA': (('5', 0.94), ('5', 0.94))}
        profile = expected
        result = strmarker.get_modal_profile(data, 'AB')
        for i in expected:
            self.assertEqual(result[i][0][0], expected[i][0][0])
            self.assertEqual(result[i][1][0], expected[i][1][0])
            self.assertAlmostEqual(result[i][0][1], expected[i][0][1])
            self.assertAlmostEqual(result[i][1][1], expected[i][1][1])
        expected = 2 * 0.1775 * 0.165 * 2 * 0.35 * 0.2525 * 2 * 0.3425 * 0.2975 * 0.94 * 0.94
        result = strmarker.calc_profile_match_probability(profile, 0.0)
        self.assertAlmostEqual(result, expected)

    def testTheta(self):
        """test theta"""
        profile_base = {'D3S1358': (('15', 113), ('17', 89)),
            'VWA': (('14', 34), ('17', 100)),
            'D16S539': (('11', 119), ('13', 70)),
            'D2S1338': (('24', 40), ('25', 48)),
            'D8S1179': (('11', 25), ('13', 139)),
            'D21S11': (('30', 105), ('31.2', 42)),
            'D18S51': (('14', 67), ('15', 57)),
            'D19S433': (('14', 131), ('15.2', 10)),
            'TH01': (('9', 55), ('9.3', 140)),
            'FGA': (('21', 71), ('21', 71))}

        profile = {}
        for i in profile_base:
            p = profile_base[i]
            profile[i] = ((p[0][0], float(p[0][1]) / 400), (p[1][0], float(p[1][1]) / 400))

        result = strmarker.calc_profile_match_probability(profile, 0.0)
        expected = 7.58e-14
        #print "result", result, expected
        self.assertAlmostEqual(result, expected, 16)

        result = strmarker.calc_profile_match_probability(profile, 0.01)
        expected = 2.44e-13
        #print "result 0.01", result, expected
        self.assertAlmostEqual(result, expected, 15)

        result = strmarker.calc_profile_match_probability(profile, 0.03)
        expected = 1.614e-12
        #print "result 0.03", result, expected
        self.assertAlmostEqual(result, expected, 14)

        # size bias correction
        for i in profile_base:
            p = profile_base[i]
            if p[0][1] == p[1][1]:
                profile[i] = ((p[0][0], float(p[0][1] + 4) / 404), (p[1][0], float(p[1][1] + 4) / 404))
            else:
                profile[i] = ((p[0][0], float(p[0][1] + 2) / 404), (p[1][0], float(p[1][1] + 2) / 404))

        result = strmarker.calc_profile_match_probability(profile, 0.0)
        expected = 1.4225e-13
        #print "result 0.0", result, expected
        self.assertAlmostEqual(result, expected, 17)

    def testRMPs(self):
        """test RMPs"""
        data = [{'name': 'AB', 'count': 400, 'marker': 'FGA', 'alleles': self.AB_Cau_FGA_alleles},
            {'name': 'AB', 'count': 400, 'marker': 'TH01', 'alleles': self.AB_Cau_TH01_alleles},
            {'name': 'AB', 'count': 400, 'marker': 'D16S539', 'alleles': self.AB_Cau_D16S539_alleles}]
        result = strmarker.calc_rmps(data, 'AB', 0, 0.0)
        expected = {'FGA': 0.036, 'TH01': 0.094, 'D16S539': 0.103}
        for i in expected:
            print i
            self.assertAlmostEqual(result[i], expected[i], 3)

    def testPoolAlleles(self):
        """test pool alleles"""
        # constructed counts (n=100)
        alleles = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 15, 'g': 30, 'h': 40}
        self.assertEqual(strmarker.pool_alleles(alleles, 0, 0), alleles)
        expected = {'other': 10, 'e': 5, 'f': 15, 'g': 30, 'h': 40}
        self.assertEqual(strmarker.pool_alleles(alleles, 5, 0), expected)

        # constructed frequencies (n=100)
        alleles = {'a': 0.01, 'b': 0.02, 'c': 0.03, 'd': 0.04, 'e': 0.05, 'f': 0.15, 'g': 0.30, 'h': 0.40}
        expected = alleles
        result = strmarker.pool_alleles(alleles, 0, 100)
        for i in expected:
            self.assertAlmostEqual(result[i], alleles[i])
        expected = {'other': 0.10, 'e': 0.05, 'f': 0.15, 'g': 0.30, 'h': 0.40}
        result = strmarker.pool_alleles(alleles, 5, 100)
        for i in expected:
            self.assertAlmostEqual(result[i], expected[i])

        # AB Cau D16S539 counts (n=200)
        alleles = {'8': 6, '9': 43, '10': 17, '11': 119, '12': 137, '13': 70, '14': 7, '15': 1}
        self.assertEqual(strmarker.pool_alleles(alleles, 0, 0), alleles)
        expected = {'9': 43, '10': 17, '11': 119, '12': 137, '13': 70, '14': 7, 'other': 7}
        self.assertEqual(strmarker.pool_alleles(alleles, 5, 0), expected)

        # AB Cau D16S539 frequencies (n=200)
        alleles = self.AB_Cau_D16S539_alleles
        expected = alleles
        result = strmarker.pool_alleles(alleles, 0, 200)
        for i in expected:
            self.assertAlmostEqual(result[i], expected[i])
        expected = self._fromPercentages({'9': 10.75, '10': 4.25, '11': 29.75, '12': 34.25, '13': 17.5, '14': 1.75,
            'other': 1.75})
        result = strmarker.pool_alleles(alleles, 5, 200 * 2)
        for i in expected:
            self.assertAlmostEqual(result[i], expected[i])


if __name__ == "__main__":
    unittest.main()
