"""
File: test_slicer_components.py

Tests components of the slicer
"""

import unittest
from analyzer import salami_slice_analyze

class SalamiSliceComponents(unittest.TestCase):
    def test_factor(self):
        """
        Tests factor
        """
        self.assertEqual(salami_slice_analyze.factor(2), [1, 2])
        self.assertEqual(salami_slice_analyze.factor(3), [1, 3])
        self.assertEqual(salami_slice_analyze.factor(10), [1, 2, 5])
        self.assertEqual(salami_slice_analyze.factor(30), [1, 2, 3, 5])
        self.assertEqual(salami_slice_analyze.factor(84), [1, 2, 2, 3, 7])
        self.assertEqual(salami_slice_analyze.factor(96), [1, 2, 2, 2, 2, 2, 3])

    def test_lcm(self):
        """
        Tests LCM
        """
        self.assertEqual(salami_slice_analyze.lcm((2, 3)), 6)
        self.assertEqual(salami_slice_analyze.lcm((2, 3, 5)), 30)
        self.assertEqual(salami_slice_analyze.lcm((2, 4)), 4)
        self.assertEqual(salami_slice_analyze.lcm((2, 3, 4)), 12)
        self.assertEqual(salami_slice_analyze.lcm((11, 3, 4)), 132)
        self.assertEqual(salami_slice_analyze.lcm((2, 3, 5, 7)), 210)

if __name__ == "__main__":
    unittest.main()