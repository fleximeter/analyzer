"""
File: test_slicer.py

Runs the slicer on music excerpts and verifies the results
"""

import unittest
import os
from pathlib import Path
from analyzer import salami_slice_analyze
import music21
from pctheory import pcset, pcseg, pset, pseg
import numpy as np
from decimal import Decimal
from fractions import Fraction

def csi(pseg: list):
    """
    Computes the chord spacing index of a pitch seg
    :param pseg: The pitch seg (note: it should not have duplicate pitches in it)
    :return: The chord spacing index (CSI)
    """
    if len(pseg) < 3:
        return np.nan
    else:
        ipseg = [p.p for p in pseg]
        return (np.average(ipseg) - ipseg[0]) / (ipseg[-1] - ipseg[0])

class BasicTests(unittest.TestCase):
    def test_simple_analysis1(self):
        """
        Tests a simple four-chord C major string quartet score
        """
        file_path = Path(__file__).parent / "data/test1.musicxml"
        analysis = salami_slice_analyze.analyze(file_path)
        self.assertEqual(len(analysis[0].slices), 4)
        self.assertEqual(analysis[0].slices[0].pseg, pseg.make_pseg12(-12, 0, 7, 16))
        self.assertEqual(analysis[0].slices[0].pcseg, pcseg.make_pcseg12(0, 0, 7, 4))
        self.assertEqual(analysis[0].slices[0].pset, pset.make_pset12(-12, 0, 7, 16))
        self.assertEqual(analysis[0].slices[0].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(analysis[0].slices[0].pitchseg, [-12, 0, 7, 16])
        self.assertEqual(analysis[0].slices[0].chord_spacing_contour, [2, 0, 1])
        self.assertAlmostEqual(analysis[0].slices[0].chord_spacing_index, csi(analysis[0].slices[0].pseg))
        
        self.assertEqual(analysis[0].slices[1].pseg, pseg.make_pseg12(-7, -3, 5, 14))
        self.assertEqual(analysis[0].slices[1].pcseg, pcseg.make_pcseg12(5, 9, 5, 2))
        self.assertEqual(analysis[0].slices[1].pset, pset.make_pset12(-7, -3, 5, 14))
        self.assertEqual(analysis[0].slices[1].pcset, pcset.make_pcset12(2, 5, 9))
        self.assertEqual(analysis[0].slices[1].pitchseg, [-7, -3, 5, 14])
        self.assertEqual(analysis[0].slices[1].chord_spacing_contour, [0, 1, 2])
        self.assertAlmostEqual(analysis[0].slices[1].chord_spacing_index, csi(analysis[0].slices[1].pseg))
        
        self.assertEqual(analysis[0].slices[2].pseg, pseg.make_pseg12(-5, 2, 11))
        self.assertEqual(analysis[0].slices[2].pcseg, pcseg.make_pcseg12(7, 2, 11))
        self.assertEqual(analysis[0].slices[2].pset, pset.make_pset12(-5, 2, 11))
        self.assertEqual(analysis[0].slices[2].pcset, pcset.make_pcset12(2, 7, 11))
        self.assertEqual(analysis[0].slices[2].pitchseg, [-5, -5, 2, 11])
        self.assertEqual(analysis[0].slices[2].chord_spacing_contour, [0, 1])
        self.assertAlmostEqual(analysis[0].slices[2].chord_spacing_index, csi(analysis[0].slices[2].pseg))

        self.assertEqual(analysis[0].slices[3].pseg, pseg.make_pseg12(-12, -5, 4, 12))
        self.assertEqual(analysis[0].slices[3].pcseg, pcseg.make_pcseg12(0, 7, 4, 0))
        self.assertEqual(analysis[0].slices[3].pset, pset.make_pset12(-12, -5, 4, 12))
        self.assertEqual(analysis[0].slices[3].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(analysis[0].slices[3].pitchseg, [-12, -5, 4, 12])
        self.assertEqual(analysis[0].slices[3].chord_spacing_contour, [0, 2, 1])
        self.assertAlmostEqual(analysis[0].slices[3].chord_spacing_index, csi(analysis[0].slices[3].pseg))

        self.assertEqual(analysis[0].slices[0].duration, 2)
        self.assertEqual(analysis[0].slices[1].duration, 2)
        self.assertEqual(analysis[0].slices[2].duration, 2)
        self.assertEqual(analysis[0].slices[3].duration, 2)

        self.assertAlmostEqual(analysis[0].chord_spacing_index_avg, np.average(
            [csi(analysis[0].slices[0].pseg), csi(analysis[0].slices[1].pseg), 
            csi(analysis[0].slices[2].pseg), csi(analysis[0].slices[3].pseg)]
        ))
    
    def test_simple_analysis2(self):
        """
        Tests a more complex version of 1 with polyrhythm
        """
        file_paths = [
            Path(__file__).parent / "data/test2.musicxml",
            Path(__file__).parent / "data/test3.musicxml",
            Path(__file__).parent / "data/test4.musicxml",
        ]
        for file_path in file_paths:
            analysis = salami_slice_analyze.analyze(file_path)
            self.assertEqual(len(analysis[0].slices), 8)
            self.assertEqual(analysis[0].slices[0].pseg, pseg.make_pseg12(-12, 0, 7, 16))
            self.assertEqual(analysis[0].slices[0].pcseg, pcseg.make_pcseg12(0, 0, 7, 4))
            self.assertEqual(analysis[0].slices[0].pset, pset.make_pset12(-12, 0, 7, 16))
            self.assertEqual(analysis[0].slices[0].pcset, pcset.make_pcset12(0, 4, 7))
            self.assertEqual(analysis[0].slices[0].pitchseg, [-12, 0, 7, 16])
            self.assertEqual(analysis[0].slices[0].chord_spacing_contour, [2, 0, 1])
            self.assertAlmostEqual(analysis[0].slices[0].chord_spacing_index, csi(analysis[0].slices[0].pseg))
            
            self.assertEqual(analysis[0].slices[1].pseg, pseg.make_pseg12(-7, -3, 5, 14))
            self.assertEqual(analysis[0].slices[1].pcseg, pcseg.make_pcseg12(5, 9, 5, 2))
            self.assertEqual(analysis[0].slices[1].pset, pset.make_pset12(-7, -3, 5, 14))
            self.assertEqual(analysis[0].slices[1].pcset, pcset.make_pcset12(2, 5, 9))
            self.assertEqual(analysis[0].slices[1].pitchseg, [-7, -3, 5, 14])
            self.assertEqual(analysis[0].slices[1].chord_spacing_contour, [0, 1, 2])
            self.assertAlmostEqual(analysis[0].slices[1].chord_spacing_index, csi(analysis[0].slices[1].pseg))
            
            self.assertEqual(analysis[0].slices[2].pseg, pseg.make_pseg12(-5, 2, 12))
            self.assertEqual(analysis[0].slices[2].pcseg, pcseg.make_pcseg12(7, 2, 0))
            self.assertEqual(analysis[0].slices[2].pset, pset.make_pset12(-5, 2, 12))
            self.assertEqual(analysis[0].slices[2].pcset, pcset.make_pcset12(2, 7, 0))
            self.assertEqual(analysis[0].slices[2].pitchseg, [-5, -5, 2, 12])
            self.assertEqual(analysis[0].slices[2].chord_spacing_contour, [0, 1])
            self.assertAlmostEqual(analysis[0].slices[2].chord_spacing_index, csi(analysis[0].slices[2].pseg))

            self.assertEqual(analysis[0].slices[3].pseg, pseg.make_pseg12(-5, 2, 11))
            self.assertEqual(analysis[0].slices[3].pcseg, pcseg.make_pcseg12(7, 2, 11))
            self.assertEqual(analysis[0].slices[3].pset, pset.make_pset12(-5, 2, 11))
            self.assertEqual(analysis[0].slices[3].pcset, pcset.make_pcset12(2, 7, 11))
            self.assertEqual(analysis[0].slices[3].pitchseg, [-5, -5, 2, 11])
            self.assertEqual(analysis[0].slices[3].chord_spacing_contour, [0, 1])
            self.assertAlmostEqual(analysis[0].slices[3].chord_spacing_index, csi(analysis[0].slices[3].pseg))

            self.assertEqual(analysis[0].slices[4].pseg, pseg.make_pseg12(-5, 2, 9))
            self.assertEqual(analysis[0].slices[4].pcseg, pcseg.make_pcseg12(7, 2, 9))
            self.assertEqual(analysis[0].slices[4].pset, pset.make_pset12(-5, 2, 9))
            self.assertEqual(analysis[0].slices[4].pcset, pcset.make_pcset12(2, 7, 9))
            self.assertEqual(analysis[0].slices[4].pitchseg, [-5, -5, 2, 9])
            self.assertEqual(analysis[0].slices[4].chord_spacing_contour, [0, 0])
            self.assertAlmostEqual(analysis[0].slices[4].chord_spacing_index, csi(analysis[0].slices[4].pseg))

            self.assertEqual(analysis[0].slices[5].pseg, pseg.make_pseg12(-5, 5, 9))
            self.assertEqual(analysis[0].slices[5].pcseg, pcseg.make_pcseg12(7, 5, 9))
            self.assertEqual(analysis[0].slices[5].pset, pset.make_pset12(-5, 5, 9))
            self.assertEqual(analysis[0].slices[5].pcset, pcset.make_pcset12(5, 7, 9))
            self.assertEqual(analysis[0].slices[5].pitchseg, [-5, -5, 5, 9])
            self.assertEqual(analysis[0].slices[5].chord_spacing_contour, [1, 0])
            self.assertAlmostEqual(analysis[0].slices[5].chord_spacing_index, csi(analysis[0].slices[5].pseg))

            self.assertEqual(analysis[0].slices[6].pseg, pseg.make_pseg12(-5, 5, 11))
            self.assertEqual(analysis[0].slices[6].pcseg, pcseg.make_pcseg12(7, 5, 11))
            self.assertEqual(analysis[0].slices[6].pset, pset.make_pset12(-5, 5, 11))
            self.assertEqual(analysis[0].slices[6].pcset, pcset.make_pcset12(5, 7, 11))
            self.assertEqual(analysis[0].slices[6].pitchseg, [-5, -5, 5, 11])
            self.assertEqual(analysis[0].slices[6].chord_spacing_contour, [1, 0])
            self.assertAlmostEqual(analysis[0].slices[6].chord_spacing_index, csi(analysis[0].slices[6].pseg))

            self.assertEqual(analysis[0].slices[7].pseg, pseg.make_pseg12(-12, -5, 4, 12))
            self.assertEqual(analysis[0].slices[7].pcseg, pcseg.make_pcseg12(0, 7, 4, 0))
            self.assertEqual(analysis[0].slices[7].pset, pset.make_pset12(-12, -5, 4, 12))
            self.assertEqual(analysis[0].slices[7].pcset, pcset.make_pcset12(0, 4, 7))
            self.assertEqual(analysis[0].slices[7].pitchseg, [-12, -5, 4, 12])
            self.assertEqual(analysis[0].slices[7].chord_spacing_contour, [0, 2, 1])
            self.assertAlmostEqual(analysis[0].slices[7].chord_spacing_index, csi(analysis[0].slices[7].pseg))

            self.assertEqual(analysis[0].slices[0].duration, 2)
            self.assertEqual(analysis[0].slices[1].duration, 2)
            self.assertEqual(analysis[0].slices[2].duration, 1)
            self.assertAlmostEqual(analysis[0].slices[3].duration, Decimal(1/3))
            self.assertAlmostEqual(analysis[0].slices[4].duration, Decimal(1/6))
            self.assertAlmostEqual(analysis[0].slices[5].duration, Decimal(1/6))
            self.assertAlmostEqual(analysis[0].slices[6].duration, Decimal(1/3))
            self.assertEqual(analysis[0].slices[7].duration, 2)

            self.assertAlmostEqual(analysis[0].chord_spacing_index_avg, 
                (csi(analysis[0].slices[0].pseg) * 2 + 
                csi(analysis[0].slices[1].pseg) * 2 + 
                csi(analysis[0].slices[2].pseg) * 1 + 
                csi(analysis[0].slices[3].pseg) * 1/3 + 
                csi(analysis[0].slices[4].pseg) * 1/6 + 
                csi(analysis[0].slices[5].pseg) * 1/6 + 
                csi(analysis[0].slices[6].pseg) * 1/3 + 
                csi(analysis[0].slices[7].pseg) * 2
                ) / 8
            )

class CarterTests(unittest.TestCase):
    def test_carter1(self):
        """
        Tests measure 310 of Carter V
        """
        results = salami_slice_analyze.analyze(Path(__file__).parent / "data/test_carter1.musicxml")
        slices = results[0].slices
        self.assertEqual(slices[0].pseg, pseg.make_pseg12(6, 8))
        self.assertEqual(slices[1].pseg, pseg.make_pseg12(3, 8, 11))
        self.assertEqual(slices[2].pseg, pseg.make_pseg12(3, 7, 8, 11))
        self.assertEqual(slices[3].pseg, pseg.make_pseg12(3, 7, 8, 9, 11, 16))
        self.assertEqual(slices[4].pseg, pseg.make_pseg12(0, 8, 9, 11, 16))
        self.assertEqual(slices[5].pseg, pseg.make_pseg12(0, 8, 11))
        self.assertEqual(slices[6].pseg, pseg.make_pseg12(0, 6, 8, 11, 13))
        self.assertEqual(slices[7].pseg, pseg.make_pseg12(0, 3, 6, 11, 13))
        self.assertEqual(slices[8].pseg, pseg.make_pseg12(0, 3, 6, 7, 11, 13))
        self.assertEqual(len(slices), 9)
        self.assertEqual(slices[0].sc_name, "(2-2)[02]")
        self.assertEqual(slices[1].sc_name, "(3-11)[037]")
        self.assertEqual(slices[2].sc_name, "(4-19)[0148]")
        self.assertEqual(slices[3].sc_name, "(6-16)[014568]")
        self.assertEqual(slices[4].sc_name, "(5-Z17)[01348]")
        self.assertEqual(slices[5].sc_name, "(3-3)[014]")
        self.assertEqual(slices[6].sc_name, "(5-14)[01257]")
        self.assertEqual(slices[7].sc_name, "(5-Z36)[01247]")
        self.assertEqual(slices[8].sc_name, "(6-Z17)[012478]")

        self.assertAlmostEqual(slices[0].duration, Decimal(0.5))
        self.assertAlmostEqual(slices[1].duration, Decimal(0.625))
        self.assertAlmostEqual(slices[2].duration, Decimal(0.208333333333))
        self.assertAlmostEqual(slices[3].duration, Decimal(0.166666666667))
        self.assertAlmostEqual(slices[4].duration, Decimal(0.166666666667))
        self.assertAlmostEqual(slices[5].duration, Decimal(0.733333333333))
        self.assertAlmostEqual(slices[6].duration, Decimal(0.1))
        self.assertAlmostEqual(slices[7].duration, Decimal(0.3))
        self.assertAlmostEqual(slices[8].duration, Decimal(1.2))
        self.assertEqual(slices[0].quarter_duration, Fraction(1, 2))
        self.assertEqual(slices[1].quarter_duration, Fraction(5, 8))
        self.assertEqual(slices[2].quarter_duration, Fraction(5, 24))
        self.assertEqual(slices[3].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[4].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[5].quarter_duration, Fraction(11, 15))
        self.assertEqual(slices[6].quarter_duration, Fraction(1, 10))
        self.assertEqual(slices[7].quarter_duration, Fraction(3, 10))
        self.assertEqual(slices[8].quarter_duration, Fraction(6, 5))

        self.assertEqual(slices[0].chord_spacing_contour, [0])
        self.assertEqual(slices[1].chord_spacing_contour, [1, 0])
        self.assertEqual(slices[2].chord_spacing_contour, [2, 0, 1])
        self.assertEqual(slices[3].chord_spacing_contour, [2, 0, 0, 1, 3])
        self.assertEqual(slices[4].chord_spacing_contour, [3, 0, 1, 2])
        self.assertEqual(slices[5].chord_spacing_contour, [1, 0])
        self.assertEqual(slices[6].chord_spacing_contour, [2, 0, 1, 0])
        self.assertEqual(slices[7].chord_spacing_contour, [1, 1, 2, 0])
        self.assertEqual(slices[8].chord_spacing_contour, [2, 2, 0, 3, 1])

        self.assertTrue(np.isnan(slices[0].chord_spacing_index))
        self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))
        self.assertAlmostEqual(slices[2].chord_spacing_index, csi(slices[2].pseg))
        self.assertAlmostEqual(slices[3].chord_spacing_index, csi(slices[3].pseg))
        self.assertAlmostEqual(slices[4].chord_spacing_index, csi(slices[4].pseg))
        self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))
        self.assertAlmostEqual(slices[6].chord_spacing_index, csi(slices[6].pseg))
        self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))
        self.assertAlmostEqual(slices[8].chord_spacing_index, csi(slices[8].pseg))

    def test_carter2(self):
        """
        Tests measure 130 of Carter V
        """
        results = salami_slice_analyze.analyze(Path(__file__).parent / "data/test_carter2.musicxml")
        slices = results[0].slices
        self.assertEqual(slices[0].pseg, pseg.make_pseg12(3, 7))
        self.assertEqual(slices[1].pseg, pseg.make_pseg12(3, 7, 13))
        self.assertEqual(slices[2].pseg, pseg.make_pseg12(13))
        self.assertEqual(slices[3].pseg, pseg.make_pseg12(12))
        self.assertEqual(slices[4].pseg, [])
        self.assertEqual(slices[5].pseg, pseg.make_pseg12(-4, 5, 14))
        self.assertEqual(slices[6].pseg, [])
        self.assertEqual(slices[7].pseg, pseg.make_pseg12(-5, 12, 15, 22))
        self.assertEqual(slices[8].pseg, pseg.make_pseg12(22))
        self.assertEqual(len(slices), 9)

        self.assertEqual(slices[0].sc_name, "(2-4)[04]")
        self.assertEqual(slices[1].sc_name, "(3-8)[026]")
        self.assertEqual(slices[2].sc_name, "(1-1)[0]")
        self.assertEqual(slices[3].sc_name, "(1-1)[0]")
        self.assertEqual(slices[4].sc_name, "(0-1)[]")
        self.assertEqual(slices[5].sc_name, "(3-10)[036]")
        self.assertEqual(slices[6].sc_name, "(0-1)[]")
        self.assertEqual(slices[7].sc_name, "(4-26)[0358]")
        self.assertEqual(slices[8].sc_name, "(1-1)[0]")

        self.assertEqual(slices[0].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[1].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[2].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[3].quarter_duration, Fraction(1, 1))
        self.assertEqual(slices[4].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[5].quarter_duration, Fraction(1, 1))
        self.assertEqual(slices[6].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[7].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[8].quarter_duration, Fraction(1, 3))
        self.assertAlmostEqual(slices[0].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[1].duration, Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(slices[2].duration, Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(slices[3].duration, Decimal(60/115.2))
        self.assertAlmostEqual(slices[4].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[5].duration, Decimal(60/115.2))
        self.assertAlmostEqual(slices[6].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[7].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[8].duration, Decimal((1/3) * (60/115.2)))

        self.assertEqual(slices[0].chord_spacing_contour, [0])
        self.assertEqual(slices[1].chord_spacing_contour, [0, 1])
        self.assertEqual(slices[2].chord_spacing_contour, [])
        self.assertEqual(slices[3].chord_spacing_contour, [])
        self.assertEqual(slices[4].chord_spacing_contour, [])
        self.assertEqual(slices[5].chord_spacing_contour, [0, 0])
        self.assertEqual(slices[6].chord_spacing_contour, [])
        self.assertEqual(slices[7].chord_spacing_contour, [2, 0, 1])
        self.assertEqual(slices[8].chord_spacing_contour, [])

        self.assertTrue(np.isnan(slices[0].chord_spacing_index))
        self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))
        self.assertTrue(np.isnan(slices[2].chord_spacing_index))
        self.assertTrue(np.isnan(slices[3].chord_spacing_index))
        self.assertTrue(np.isnan(slices[4].chord_spacing_index))
        self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))
        self.assertTrue(np.isnan(slices[6].chord_spacing_index))
        self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))
        self.assertTrue(np.isnan(slices[8].chord_spacing_index))

if __name__ == "__main__":
    unittest.main()