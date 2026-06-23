import unittest
import os
from pathlib import Path
from analyzer import salami_slice_analyze
import music21
from pctheory import pcset, pcseg, pset, pseg

class BasicTest1(unittest.TestCase):
    def test_simple_analysis(self):
        file_path = Path(__file__).parent / "data/test1.musicxml"
        analysis = salami_slice_analyze.analyze(file_path)
        self.assertEqual(len(analysis[0].slices), 4)
        self.assertEqual(analysis[0].slices[0].pseg, pseg.make_pseg12(-12, 0, 7, 16))
        self.assertEqual(analysis[0].slices[0].pcseg, pcseg.make_pcseg12(0, 0, 7, 4))
        self.assertEqual(analysis[0].slices[0].pset, pset.make_pset12(-12, 0, 7, 16))
        self.assertEqual(analysis[0].slices[0].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(analysis[0].slices[0].pitchseg, [-12, 0, 7, 16])
        self.assertEqual(analysis[0].slices[0].chord_spacing_contour, [2, 0, 1])
        
        self.assertEqual(analysis[0].slices[1].pseg, pseg.make_pseg12(-7, -3, 5, 14))
        self.assertEqual(analysis[0].slices[1].pcseg, pcseg.make_pcseg12(5, 9, 5, 2))
        self.assertEqual(analysis[0].slices[1].pset, pset.make_pset12(-7, -3, 5, 14))
        self.assertEqual(analysis[0].slices[1].pcset, pcset.make_pcset12(2, 5, 9))
        self.assertEqual(analysis[0].slices[1].pitchseg, [-7, -3, 5, 14])
        self.assertEqual(analysis[0].slices[1].chord_spacing_contour, [0, 1, 2])
        
        self.assertEqual(analysis[0].slices[2].pseg, pseg.make_pseg12(-5, 2, 11))
        self.assertEqual(analysis[0].slices[2].pcseg, pcseg.make_pcseg12(7, 2, 11))
        self.assertEqual(analysis[0].slices[2].pset, pset.make_pset12(-5, 2, 11))
        self.assertEqual(analysis[0].slices[2].pcset, pcset.make_pcset12(2, 7, 11))
        self.assertEqual(analysis[0].slices[2].pitchseg, [-5, -5, 2, 11])
        self.assertEqual(analysis[0].slices[2].chord_spacing_contour, [0, 1])
        
        self.assertEqual(analysis[0].slices[3].pseg, pseg.make_pseg12(-12, -5, 4, 12))
        self.assertEqual(analysis[0].slices[3].pcseg, pcseg.make_pcseg12(0, 7, 4, 0))
        self.assertEqual(analysis[0].slices[3].pset, pset.make_pset12(-12, -5, 4, 12))
        self.assertEqual(analysis[0].slices[3].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(analysis[0].slices[3].pitchseg, [-12, -5, 4, 12])
        self.assertEqual(analysis[0].slices[3].chord_spacing_contour, [0, 2, 1])
        
        self.assertEqual(analysis[0].slices[0].duration, 2)
        self.assertEqual(analysis[0].slices[1].duration, 2)
        self.assertEqual(analysis[0].slices[2].duration, 2)
        self.assertEqual(analysis[0].slices[3].duration, 2)

if __name__ == "__main__":
    unittest.main()