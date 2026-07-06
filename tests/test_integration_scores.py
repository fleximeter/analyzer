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
from pctheory.pitch import Pitch, PitchClass
from pctheory.pcset import SetClass

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
        analysis = salami_slice_analyze.analyze(file_path)[0]
        slices = analysis.slices

        self.assertEqual(len(slices), 4)
        self.assertEqual(slices[0].pseg, pseg.make_pseg12(-12, 0, 7, 16))
        self.assertEqual(slices[0].pcseg, pcseg.make_pcseg12(0, 0, 7, 4))
        self.assertEqual(slices[0].pset, pset.make_pset12(-12, 0, 7, 16))
        self.assertEqual(slices[0].psets, [{Pitch(16),}, {Pitch(7),}, {Pitch(0),}, {Pitch(-12),}])
        self.assertEqual(slices[0].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(slices[0].pitchseg, [-12, 0, 7, 16])
        self.assertEqual(slices[0].chord_spacing_contour, [2, 0, 1])
        self.assertAlmostEqual(slices[0].chord_spacing_index, csi(slices[0].pseg))
        
        self.assertEqual(slices[1].pseg, pseg.make_pseg12(-7, -3, 5, 14))
        self.assertEqual(slices[1].pcseg, pcseg.make_pcseg12(5, 9, 5, 2))
        self.assertEqual(slices[1].pset, pset.make_pset12(-7, -3, 5, 14))
        self.assertEqual(slices[1].psets, [{Pitch(14),}, {Pitch(5),}, {Pitch(-3),}, {Pitch(-7),}])
        self.assertEqual(slices[1].pcset, pcset.make_pcset12(2, 5, 9))
        self.assertEqual(slices[1].pitchseg, [-7, -3, 5, 14])
        self.assertEqual(slices[1].chord_spacing_contour, [0, 1, 2])
        self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))
        
        self.assertEqual(slices[2].pseg, pseg.make_pseg12(-5, 2, 11))
        self.assertEqual(slices[2].pcseg, pcseg.make_pcseg12(7, 2, 11))
        self.assertEqual(slices[2].pset, pset.make_pset12(-5, 2, 11))
        self.assertEqual(slices[2].psets, [{Pitch(11),}, {Pitch(2),}, {Pitch(-5),}, {Pitch(-5),}])
        self.assertEqual(slices[2].pcset, pcset.make_pcset12(2, 7, 11))
        self.assertEqual(slices[2].pitchseg, [-5, -5, 2, 11])
        self.assertEqual(slices[2].chord_spacing_contour, [0, 1])
        self.assertAlmostEqual(slices[2].chord_spacing_index, csi(slices[2].pseg))

        self.assertEqual(slices[3].pseg, pseg.make_pseg12(-12, -5, 4, 12))
        self.assertEqual(slices[3].pcseg, pcseg.make_pcseg12(0, 7, 4, 0))
        self.assertEqual(slices[3].pset, pset.make_pset12(-12, -5, 4, 12))
        self.assertEqual(slices[3].psets, [{Pitch(12),}, {Pitch(4),}, {Pitch(-5),}, {Pitch(-12),}])
        self.assertEqual(slices[3].pcset, pcset.make_pcset12(0, 4, 7))
        self.assertEqual(slices[3].pitchseg, [-12, -5, 4, 12])
        self.assertEqual(slices[3].chord_spacing_contour, [0, 2, 1])
        self.assertAlmostEqual(slices[3].chord_spacing_index, csi(slices[3].pseg))

        self.assertEqual(slices[0].duration, 2)
        self.assertEqual(slices[1].duration, 2)
        self.assertEqual(slices[2].duration, 2)
        self.assertEqual(slices[3].duration, 2)

        self.assertAlmostEqual(analysis.chord_spacing_index_avg, np.average(
            [csi(slices[0].pseg), csi(slices[1].pseg), 
            csi(slices[2].pseg), csi(slices[3].pseg)]
        ))

        self.assertEqual(analysis.duration_avg, 2)
        self.assertEqual(analysis.duration, 8)
        self.assertEqual(analysis.quarter_duration, 8)
        self.assertEqual(analysis.measure_num_first, 1)
        self.assertEqual(analysis.measure_num_last, 2)
        self.assertEqual(analysis.num_measures, 2)
        self.assertEqual(analysis.num_voices, 4)
        self.assertEqual(analysis.max_pitch_count_with_duplicates, 4)
        self.assertEqual(analysis.lower_bound, -12)
        self.assertEqual(analysis.upper_bound, 16)
        self.assertEqual(analysis.start_time, 0)
        self.assertEqual(analysis.pc_duration, {
            0: Decimal(4), 2: Decimal(4), 
            4: Decimal(4), 5: Decimal(2),
            7: Decimal(6), 9: Decimal(2), 11: Decimal(2)
        })
        self.assertEqual(analysis.pitch_duration, {
            -12: Decimal(4), -7: Decimal(2), -5: Decimal(4),
            -3: Decimal(2), 0: Decimal(2), 2: Decimal(2),
            4: Decimal(2), 5: Decimal(2), 7: Decimal(2),
            11: Decimal(2), 12: Decimal(2), 14: Decimal(2), 16: Decimal(2)
        })
        # Right now this is a little weird. Basically if a pc reoccurs in
        # the very next slice, it isn't counted as a new occurrence for
        # frequency. This holds even if it's a different pitch in a different
        # octave. It's not clear if this behavior should be this way or not,
        # but for now it's probably best not to rely on this statistic in a paper.
        self.assertEqual(analysis.pc_frequency, {
            0: 2, 2: 1, 4: 2, 5: 1, 7: 2, 9: 1, 11: 1
        })
        # Same holds here. If the pitch occurs in the next slice (even if it's in
        # a different instrument), it isn't counted again.
        self.assertEqual(analysis.pitch_frequency, {
            -12: 2, -7: 1, -5: 1, -3: 1, 0: 1, 2: 1, 4: 1,
            5: 1, 7: 1, 11: 1, 12: 1, 14: 1, 16: 1
        })

        self.assertEqual(analysis.pcsc_frequency, {SetClass("[037]").name_morris: 4})
        self.assertEqual(analysis.pcsc_duration, {SetClass("[037]").name_morris: 8})
        self.assertEqual(analysis.pset_frequency, {
            "{-12, 0, 7, 16}": 1,
            "{-7, -3, 5, 14}": 1,
            "{-5, 2, 11}": 1,
            "{-12, -5, 4, 12}": 1
        })
        self.assertEqual(analysis.pset_duration, {
            "{-12, 0, 7, 16}": 2,
            "{-7, -3, 5, 14}": 2,
            "{-5, 2, 11}": 2,
            "{-12, -5, 4, 12}": 2
        })
        self.assertEqual(analysis.psc_frequency, {
            "[12, 7, 9]": 1,
            "[4, 8, 9]": 1,
            "[7, 9]": 1,
            "[7, 9, 8]": 1
        })
        self.assertEqual(analysis.psc_duration, {
            "[12, 7, 9]": 2,
            "[4, 8, 9]": 2,
            "[7, 9]": 2,
            "[7, 9, 8]": 2
        })
        self.assertEqual(analysis.chord_spacing_contour_frequency, {
            "<2, 0, 1>": 1,
            "<0, 1, 2>": 1,
            "<0, 1>": 1,
            "<0, 2, 1>": 1
        })
        self.assertEqual(analysis.chord_spacing_contour_duration, {
            "<2, 0, 1>": 2,
            "<0, 1, 2>": 2,
            "<0, 1>": 2,
            "<0, 2, 1>": 2
        })

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
            analysis = salami_slice_analyze.analyze(file_path)[0]
            slices = analysis.slices
            self.assertEqual(len(slices), 8)
            self.assertEqual(slices[0].pseg, pseg.make_pseg12(-12, 0, 7, 16))
            self.assertEqual(slices[0].pcseg, pcseg.make_pcseg12(0, 0, 7, 4))
            self.assertEqual(slices[0].pset, pset.make_pset12(-12, 0, 7, 16))            
            self.assertEqual(slices[0].pcset, pcset.make_pcset12(0, 4, 7))
            self.assertEqual(slices[0].pitchseg, [-12, 0, 7, 16])
            self.assertEqual(slices[0].chord_spacing_contour, [2, 0, 1])
            self.assertAlmostEqual(slices[0].chord_spacing_index, csi(slices[0].pseg))
            
            self.assertEqual(slices[1].pseg, pseg.make_pseg12(-7, -3, 5, 14))
            self.assertEqual(slices[1].pcseg, pcseg.make_pcseg12(5, 9, 5, 2))
            self.assertEqual(slices[1].pset, pset.make_pset12(-7, -3, 5, 14))
            self.assertEqual(slices[1].pcset, pcset.make_pcset12(2, 5, 9))
            self.assertEqual(slices[1].pitchseg, [-7, -3, 5, 14])
            self.assertEqual(slices[1].chord_spacing_contour, [0, 1, 2])
            self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))
            
            self.assertEqual(slices[2].pseg, pseg.make_pseg12(-5, 2, 12))
            self.assertEqual(slices[2].pcseg, pcseg.make_pcseg12(7, 2, 0))
            self.assertEqual(slices[2].pset, pset.make_pset12(-5, 2, 12))
            self.assertEqual(slices[2].pcset, pcset.make_pcset12(2, 7, 0))
            self.assertEqual(slices[2].pitchseg, [-5, -5, 2, 12])
            self.assertEqual(slices[2].chord_spacing_contour, [0, 1])
            self.assertAlmostEqual(slices[2].chord_spacing_index, csi(slices[2].pseg))

            self.assertEqual(slices[3].pseg, pseg.make_pseg12(-5, 2, 11))
            self.assertEqual(slices[3].pcseg, pcseg.make_pcseg12(7, 2, 11))
            self.assertEqual(slices[3].pset, pset.make_pset12(-5, 2, 11))
            self.assertEqual(slices[3].pcset, pcset.make_pcset12(2, 7, 11))
            self.assertEqual(slices[3].pitchseg, [-5, -5, 2, 11])
            self.assertEqual(slices[3].chord_spacing_contour, [0, 1])
            self.assertAlmostEqual(slices[3].chord_spacing_index, csi(slices[3].pseg))

            self.assertEqual(slices[4].pseg, pseg.make_pseg12(-5, 2, 9))
            self.assertEqual(slices[4].pcseg, pcseg.make_pcseg12(7, 2, 9))
            self.assertEqual(slices[4].pset, pset.make_pset12(-5, 2, 9))
            self.assertEqual(slices[4].pcset, pcset.make_pcset12(2, 7, 9))
            self.assertEqual(slices[4].pitchseg, [-5, -5, 2, 9])
            self.assertEqual(slices[4].chord_spacing_contour, [0, 0])
            self.assertAlmostEqual(slices[4].chord_spacing_index, csi(slices[4].pseg))

            self.assertEqual(slices[5].pseg, pseg.make_pseg12(-5, 5, 9))
            self.assertEqual(slices[5].pcseg, pcseg.make_pcseg12(7, 5, 9))
            self.assertEqual(slices[5].pset, pset.make_pset12(-5, 5, 9))
            self.assertEqual(slices[5].pcset, pcset.make_pcset12(5, 7, 9))
            self.assertEqual(slices[5].pitchseg, [-5, -5, 5, 9])
            self.assertEqual(slices[5].chord_spacing_contour, [1, 0])
            self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))

            self.assertEqual(slices[6].pseg, pseg.make_pseg12(-5, 5, 11))
            self.assertEqual(slices[6].pcseg, pcseg.make_pcseg12(7, 5, 11))
            self.assertEqual(slices[6].pset, pset.make_pset12(-5, 5, 11))
            self.assertEqual(slices[6].pcset, pcset.make_pcset12(5, 7, 11))
            self.assertEqual(slices[6].pitchseg, [-5, -5, 5, 11])
            self.assertEqual(slices[6].chord_spacing_contour, [1, 0])
            self.assertAlmostEqual(slices[6].chord_spacing_index, csi(slices[6].pseg))

            self.assertEqual(slices[7].pseg, pseg.make_pseg12(-12, -5, 4, 12))
            self.assertEqual(slices[7].pcseg, pcseg.make_pcseg12(0, 7, 4, 0))
            self.assertEqual(slices[7].pset, pset.make_pset12(-12, -5, 4, 12))
            self.assertEqual(slices[7].pcset, pcset.make_pcset12(0, 4, 7))
            self.assertEqual(slices[7].pitchseg, [-12, -5, 4, 12])
            self.assertEqual(slices[7].chord_spacing_contour, [0, 2, 1])
            self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))

            self.assertEqual(slices[0].duration, 2)
            self.assertEqual(slices[1].duration, 2)
            self.assertEqual(slices[2].duration, 1)
            self.assertAlmostEqual(slices[3].duration, Decimal(1/3))
            self.assertAlmostEqual(slices[4].duration, Decimal(1/6))
            self.assertAlmostEqual(slices[5].duration, Decimal(1/6))
            self.assertAlmostEqual(slices[6].duration, Decimal(1/3))
            self.assertEqual(slices[7].duration, 2)

            self.assertAlmostEqual(analysis.chord_spacing_index_avg, 
                (csi(slices[0].pseg) * 2 + 
                csi(slices[1].pseg) * 2 + 
                csi(slices[2].pseg) * 1 + 
                csi(slices[3].pseg) * 1/3 + 
                csi(slices[4].pseg) * 1/6 + 
                csi(slices[5].pseg) * 1/6 + 
                csi(slices[6].pseg) * 1/3 + 
                csi(slices[7].pseg) * 2
                ) / 8
            )

            self.assertAlmostEqual(analysis.duration_avg, np.average([2,2,1,1/3,1/6,1/6,1/3,2]))
            self.assertEqual(analysis.duration, 8)
            self.assertEqual(analysis.quarter_duration, 8)
            self.assertEqual(analysis.measure_num_first, 1)
            self.assertEqual(analysis.measure_num_last, 2)
            self.assertEqual(analysis.num_measures, 2)
            self.assertEqual(analysis.max_pitch_count_with_duplicates, 4)
            self.assertEqual(analysis.lower_bound, -12)
            self.assertEqual(analysis.upper_bound, 16)
            self.assertEqual(analysis.start_time, 0)
            self.assertEqual(len(analysis.pc_duration), 7)
            self.assertAlmostEqual(analysis.pc_duration[0], Decimal(5))
            self.assertAlmostEqual(analysis.pc_duration[2], Decimal(3.5))
            self.assertAlmostEqual(analysis.pc_duration[4], Decimal(4))
            self.assertAlmostEqual(analysis.pc_duration[5], Decimal(2.5))
            self.assertAlmostEqual(analysis.pc_duration[7], Decimal(6))
            self.assertAlmostEqual(analysis.pc_duration[9], Decimal(7/3))
            self.assertAlmostEqual(analysis.pc_duration[11], Decimal(2/3))
            self.assertEqual(len(analysis.pitch_duration), 14)
            self.assertAlmostEqual(analysis.pitch_duration[-12], Decimal(4))
            self.assertAlmostEqual(analysis.pitch_duration[-7], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[-5], Decimal(4))
            self.assertAlmostEqual(analysis.pitch_duration[-3], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[0], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[2], Decimal(1.5))
            self.assertAlmostEqual(analysis.pitch_duration[4], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[5], Decimal(2.5))
            self.assertAlmostEqual(analysis.pitch_duration[7], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[9], Decimal(1/3))
            self.assertAlmostEqual(analysis.pitch_duration[11], Decimal(2/3))
            self.assertAlmostEqual(analysis.pitch_duration[12], Decimal(3))
            self.assertAlmostEqual(analysis.pitch_duration[14], Decimal(2))
            self.assertAlmostEqual(analysis.pitch_duration[16], Decimal(2))
            self.assertEqual(analysis.pc_frequency, {
                0: 3, 2: 1, 4: 2, 5: 2, 7: 2, 9: 2, 11: 2
            })
            self.assertEqual(analysis.pitch_frequency, {
                -12: 2, -7: 1, -5: 1, -3: 1, 0: 1, 2: 1, 4: 1,
                5: 2, 7: 1, 9: 1, 11: 2, 12: 2, 14: 1, 16: 1
            })
            self.assertEqual(analysis.pcsc_frequency, {
                SetClass("[037]").name_morris: 4,
                SetClass("[027]").name_morris: 2,
                SetClass("[024]").name_morris: 1,
                SetClass("[026]").name_morris: 1
            })
            self.assertEqual(len(analysis.pcsc_duration), 4)
            self.assertAlmostEqual(analysis.pcsc_duration[SetClass("[037]").name_morris], Decimal(19/3))
            self.assertAlmostEqual(analysis.pcsc_duration[SetClass("[027]").name_morris], Decimal(7/6))
            self.assertAlmostEqual(analysis.pcsc_duration[SetClass("[024]").name_morris], Decimal(1/6))
            self.assertAlmostEqual(analysis.pcsc_duration[SetClass("[026]").name_morris], Decimal(1/3))
            self.assertEqual(analysis.pset_frequency, {
                "{-12, 0, 7, 16}": 1,
                "{-7, -3, 5, 14}": 1,
                "{-5, 2, 12}": 1,
                "{-5, 2, 11}": 1,
                "{-5, 2, 9}": 1,
                "{-5, 5, 9}": 1,
                "{-5, 5, 11}": 1,
                "{-12, -5, 4, 12}": 1
            })
            self.assertEqual(len(analysis.pset_duration), 8)
            self.assertAlmostEqual(analysis.pset_duration["{-12, 0, 7, 16}"], Decimal(2))
            self.assertAlmostEqual(analysis.pset_duration["{-7, -3, 5, 14}"], Decimal(2))
            self.assertAlmostEqual(analysis.pset_duration["{-5, 2, 12}"], Decimal(1))
            self.assertAlmostEqual(analysis.pset_duration["{-5, 2, 11}"], Decimal(1/3))
            self.assertAlmostEqual(analysis.pset_duration["{-5, 2, 9}"], Decimal(1/6))
            self.assertAlmostEqual(analysis.pset_duration["{-5, 5, 9}"], Decimal(1/6))
            self.assertAlmostEqual(analysis.pset_duration["{-5, 5, 11}"], Decimal(1/3))
            self.assertAlmostEqual(analysis.pset_duration["{-12, -5, 4, 12}"], Decimal(2))
            self.assertEqual(analysis.psc_frequency, {
                "[12, 7, 9]": 1,
                "[4, 8, 9]": 1,
                "[7, 10]": 1,
                "[7, 9]": 1,
                "[7, 7]": 1,
                "[10, 4]": 1,
                "[10, 6]": 1,
                "[7, 9, 8]": 1
            })
            self.assertEqual(len(analysis.psc_duration), 8)
            self.assertAlmostEqual(analysis.psc_duration["[12, 7, 9]"], Decimal(2))
            self.assertAlmostEqual(analysis.psc_duration["[4, 8, 9]"], Decimal(2))
            self.assertAlmostEqual(analysis.psc_duration["[7, 10]"], Decimal(1))
            self.assertAlmostEqual(analysis.psc_duration["[7, 9]"], Decimal(1/3))
            self.assertAlmostEqual(analysis.psc_duration["[7, 7]"], Decimal(1/6))
            self.assertAlmostEqual(analysis.psc_duration["[10, 4]"], Decimal(1/6))
            self.assertAlmostEqual(analysis.psc_duration["[10, 6]"], Decimal(1/3))
            self.assertAlmostEqual(analysis.psc_duration["[7, 9, 8]"], Decimal(2))
            self.assertEqual(analysis.chord_spacing_contour_frequency, {
                "<2, 0, 1>": 1,
                "<0, 1, 2>": 1,
                "<0, 1>": 2,
                "<0, 0>": 1,
                "<1, 0>": 2,
                "<0, 2, 1>": 1
            })
            self.assertEqual(len(analysis.chord_spacing_contour_duration), 6)
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<2, 0, 1>"], Decimal(2))
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<0, 1, 2>"], Decimal(2))
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<0, 1>"], Decimal(4/3))
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<0, 0>"], Decimal(1/6))
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<1, 0>"], Decimal(0.5))
            self.assertAlmostEqual(analysis.chord_spacing_contour_duration["<0, 2, 1>"], Decimal(2))

class CarterTests(unittest.TestCase):
    def test_carter1(self):
        """
        Tests measure 310 of Carter V
        """
        results = salami_slice_analyze.analyze(Path(__file__).parent / "data/test_carter1.musicxml")[0]
        slices = results.slices
        self.assertEqual(len(slices), 9)
        
        self.assertEqual(slices[0].pseg, pseg.make_pseg12(6, 8))
        self.assertEqual(slices[0].pset, pset.make_pset12(6, 8))
        self.assertEqual(slices[0].psets, [{Pitch(6), Pitch(8)}, set(), set(), set()])
        self.assertEqual(slices[0].pcset, pcset.make_pcset12(6, 8))
        self.assertEqual(slices[0].pcseg, pcseg.make_pcseg12(6, 8))
        self.assertEqual(slices[0].pitchseg, [6, 8])
        self.assertEqual(slices[0].sc_name, "(2-2)[02]")
        self.assertEqual(slices[0].chord_spacing_contour, [0])
        self.assertTrue(np.isnan(slices[0].chord_spacing_index))

        self.assertEqual(slices[1].pseg, pseg.make_pseg12(3, 8, 11))
        self.assertEqual(slices[1].pset, pset.make_pset12(3, 8, 11))
        self.assertEqual(slices[1].psets, [pset.make_pset12(3, 8, 11), set(), set(), set()])
        self.assertEqual(slices[1].pcset, pcset.make_pcset12(3, 8, 11))
        self.assertEqual(slices[1].pcseg, pcseg.make_pcseg12(3, 8, 11))
        self.assertEqual(slices[1].pitchseg, [3, 8, 11])
        self.assertEqual(slices[1].sc_name, "(3-11)[037]")
        self.assertEqual(slices[1].chord_spacing_contour, [1, 0])
        self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))

        self.assertEqual(slices[2].pseg, pseg.make_pseg12(3, 7, 8, 11))
        self.assertEqual(slices[2].pset, pset.make_pset12(3, 7, 8, 11))
        self.assertEqual(slices[2].psets, 
            [pset.make_pset12(3, 8, 11), set(), pset.make_pset12(3, 7), set()])
        self.assertEqual(slices[2].pcset, pcset.make_pcset12(3, 7, 8, 11))
        self.assertEqual(slices[2].pcseg, pcseg.make_pcseg12(3, 7, 8, 11))
        self.assertEqual(slices[2].pitchseg, [3, 3, 7, 8, 11])
        self.assertEqual(slices[2].sc_name, "(4-19)[0148]")
        self.assertEqual(slices[2].chord_spacing_contour, [2, 0, 1])
        self.assertAlmostEqual(slices[2].chord_spacing_index, csi(slices[2].pseg))

        self.assertEqual(slices[3].pseg, pseg.make_pseg12(3, 7, 8, 9, 11, 16))
        self.assertEqual(slices[3].pset, pset.make_pset12(3, 7, 8, 9, 11, 16))
        self.assertEqual(slices[3].psets, 
            [pset.make_pset12(3, 8, 11), set(), 
            pset.make_pset12(3, 7), pset.make_pset12(9, 16)])
        self.assertEqual(slices[3].pcset, pcset.make_pcset12(3, 4, 7, 8, 9, 11))
        self.assertEqual(slices[3].pcseg, pcseg.make_pcseg12(3, 7, 8, 9, 11, 4))
        self.assertEqual(slices[3].pitchseg, [3, 3, 7, 8, 9, 11, 16])
        self.assertEqual(slices[3].sc_name, "(6-16)[014568]")
        self.assertEqual(slices[3].chord_spacing_contour, [2, 0, 0, 1, 3])
        self.assertAlmostEqual(slices[3].chord_spacing_index, csi(slices[3].pseg))

        self.assertEqual(slices[4].pseg, pseg.make_pseg12(0, 8, 9, 11, 16))
        self.assertEqual(slices[4].pset, pset.make_pset12(0, 8, 9, 11, 16))
        self.assertEqual(slices[4].psets, 
            [pset.make_pset12(0, 8, 11), set(), set(), pset.make_pset12(9, 16)])
        self.assertEqual(slices[4].pcset, pcset.make_pcset12(0, 4, 8, 9, 11))
        self.assertEqual(slices[4].pcseg, pcseg.make_pcseg12(0, 8, 9, 11, 4))
        self.assertEqual(slices[4].pitchseg, [0, 8, 9, 11, 16])
        self.assertEqual(slices[4].sc_name, "(5-Z17)[01348]")
        self.assertEqual(slices[4].chord_spacing_contour, [3, 0, 1, 2])
        self.assertAlmostEqual(slices[4].chord_spacing_index, csi(slices[4].pseg))

        self.assertEqual(slices[5].pseg, pseg.make_pseg12(0, 8, 11))
        self.assertEqual(slices[5].pset, pset.make_pset12(0, 8, 11))
        self.assertEqual(slices[5].psets, 
            [pset.make_pset12(0, 8, 11), set(), set(), set()])
        self.assertEqual(slices[5].pcset, pcset.make_pcset12(0, 8, 11))
        self.assertEqual(slices[5].pcseg, pcseg.make_pcseg12(0, 8, 11))
        self.assertEqual(slices[5].pitchseg, [0, 8, 11])
        self.assertEqual(slices[5].sc_name, "(3-3)[014]")
        self.assertEqual(slices[5].chord_spacing_contour, [1, 0])
        self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))

        self.assertEqual(slices[6].pseg, pseg.make_pseg12(0, 6, 8, 11, 13))
        self.assertEqual(slices[6].pset, pset.make_pset12(0, 6, 8, 11, 13))
        self.assertEqual(slices[6].psets, 
            [pset.make_pset12(0, 8, 11), pset.make_pset12(0, 6, 13), set(), set()])
        self.assertEqual(slices[6].pcset, pcset.make_pcset12(0, 6, 8, 11, 13))
        self.assertEqual(slices[6].pcseg, pcseg.make_pcseg12(0, 6, 8, 11, 13))
        self.assertEqual(slices[6].pitchseg, [0, 0, 6, 8, 11, 13])
        self.assertEqual(slices[6].sc_name, "(5-14)[01257]")
        self.assertEqual(slices[6].chord_spacing_contour, [2, 0, 1, 0])
        self.assertAlmostEqual(slices[6].chord_spacing_index, csi(slices[6].pseg))

        self.assertEqual(slices[7].pseg, pseg.make_pseg12(0, 3, 6, 11, 13))
        self.assertEqual(slices[7].pset, pset.make_pset12(0, 3, 6, 11, 13))
        self.assertEqual(slices[7].psets, 
            [pset.make_pset12(3, 6, 11), pset.make_pset12(0, 6, 13), set(), set()])
        self.assertEqual(slices[7].pcset, pcset.make_pcset12(0, 3, 6, 11, 13))
        self.assertEqual(slices[7].pcseg, pcseg.make_pcseg12(0, 3, 6, 11, 13))
        self.assertEqual(slices[7].pitchseg, [0, 3, 6, 6, 11, 13])
        self.assertEqual(slices[7].sc_name, "(5-Z36)[01247]")
        self.assertEqual(slices[7].chord_spacing_contour, [1, 1, 2, 0])
        self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))

        self.assertEqual(slices[8].pseg, pseg.make_pseg12(0, 3, 6, 7, 11, 13))
        self.assertEqual(slices[8].pset, pset.make_pset12(0, 3, 6, 7, 11, 13))
        self.assertEqual(slices[8].psets, 
            [pset.make_pset12(3, 6, 11), pset.make_pset12(0, 7, 13), set(), set()])
        self.assertEqual(slices[8].pcset, pcset.make_pcset12(0, 3, 6, 7, 11, 13))
        self.assertEqual(slices[8].pcseg, pcseg.make_pcseg12(0, 3, 6, 7, 11, 13))
        self.assertEqual(slices[8].pitchseg, [0, 3, 6, 7, 11, 13])
        self.assertEqual(slices[8].sc_name, "(6-Z17)[012478]")
        self.assertEqual(slices[8].chord_spacing_contour, [2, 2, 0, 3, 1])
        self.assertAlmostEqual(slices[8].chord_spacing_index, csi(slices[8].pseg))

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

        # first slice has CSI of nan
        self.assertAlmostEqual(results.chord_spacing_index_avg, (
            float(slices[1].duration) * csi(slices[1].pseg) + 
            float(slices[2].duration) * csi(slices[2].pseg) + 
            float(slices[3].duration) * csi(slices[3].pseg) + 
            float(slices[4].duration) * csi(slices[4].pseg) + 
            float(slices[5].duration) * csi(slices[5].pseg) + 
            float(slices[6].duration) * csi(slices[6].pseg) + 
            float(slices[7].duration) * csi(slices[7].pseg) + 
            float(slices[8].duration) * csi(slices[8].pseg)
        ) / float(results.duration - slices[0].duration))

        self.assertAlmostEqual(float(results.duration_avg), np.average(
            [float(s.duration) for s in slices]
        ))
        self.assertAlmostEqual(results.duration, 4.0)
        self.assertEqual(results.quarter_duration, 4)
        self.assertEqual(results.measure_num_first, 1)
        self.assertEqual(results.measure_num_last, 1)
        self.assertEqual(results.num_measures, 1)
        self.assertEqual(results.num_voices, 4)
        self.assertEqual(results.max_pitch_count_with_duplicates, 7)
        self.assertEqual(results.lower_bound, 0)
        self.assertEqual(results.upper_bound, 16)
        self.assertEqual(results.start_time, 0)
        self.assertEqual(len(results.pc_duration), 9)
        self.assertAlmostEqual(results.pc_duration[0], Decimal(2.5))
        self.assertAlmostEqual(results.pc_duration[1], Decimal(8/5))
        self.assertAlmostEqual(results.pc_duration[3], Decimal(2.5))
        self.assertAlmostEqual(results.pc_duration[4], Decimal(1/3))
        self.assertAlmostEqual(results.pc_duration[6], Decimal(8/5+0.5))
        self.assertAlmostEqual(results.pc_duration[7], Decimal(6/5+0.375))
        self.assertAlmostEqual(results.pc_duration[8], Decimal(2.5))
        self.assertAlmostEqual(results.pc_duration[9], Decimal(1/3))
        self.assertAlmostEqual(results.pc_duration[11], Decimal(3.5))
        self.assertEqual(len(results.pitch_duration), 9)
        self.assertAlmostEqual(results.pitch_duration[0], Decimal(2.5))
        self.assertAlmostEqual(results.pitch_duration[3], Decimal(2.5))
        self.assertAlmostEqual(results.pitch_duration[6], Decimal(8/5+0.5))
        self.assertAlmostEqual(results.pitch_duration[7], Decimal(6/5+0.375))
        self.assertAlmostEqual(results.pitch_duration[8], Decimal(2.5))
        self.assertAlmostEqual(results.pitch_duration[9], Decimal(1/3))
        self.assertAlmostEqual(results.pitch_duration[11], Decimal(3.5))
        self.assertAlmostEqual(results.pitch_duration[13], Decimal(8/5))
        self.assertAlmostEqual(results.pitch_duration[16], Decimal(1/3))
        self.assertEqual(results.pc_frequency, {
            0: 1, 1: 1, 3: 2, 4: 1, 6: 2, 7: 2, 8: 1, 9: 1, 11: 1
        })
        self.assertEqual(results.pitch_frequency, {
            0: 1, 3: 2, 6: 2, 7: 2, 8: 1, 9: 1, 11: 1, 13: 1, 16: 1
        })
        self.assertEqual(results.pcsc_frequency, {
            "(2-2)[02]": 1,
            "(3-11)[037]": 1,
            "(4-19)[0148]": 1,
            "(6-16)[014568]": 1,
            "(5-Z17)[01348]": 1,
            "(3-3)[014]": 1,
            "(5-14)[01257]": 1,
            "(5-Z36)[01247]": 1,
            "(6-Z17)[012478]": 1
        })
        self.assertEqual(len(results.pcsc_duration), 9)
        self.assertAlmostEqual(results.pcsc_duration["(2-2)[02]"], Decimal(0.5))
        self.assertAlmostEqual(results.pcsc_duration["(3-11)[037]"], Decimal(0.625))
        self.assertAlmostEqual(results.pcsc_duration["(4-19)[0148]"], Decimal(1/3-0.125))
        self.assertAlmostEqual(results.pcsc_duration["(6-16)[014568]"], Decimal(1/6))
        self.assertAlmostEqual(results.pcsc_duration["(5-Z17)[01348]"], Decimal(1/6))
        self.assertAlmostEqual(results.pcsc_duration["(3-3)[014]"], Decimal(1/3+0.4))
        self.assertAlmostEqual(results.pcsc_duration["(5-14)[01257]"], Decimal(0.1))
        self.assertAlmostEqual(results.pcsc_duration["(5-Z36)[01247]"], Decimal(0.3))
        self.assertAlmostEqual(results.pcsc_duration["(6-Z17)[012478]"], Decimal(6/5))
        self.assertEqual(results.psc_frequency, {
            "[2]": 1,
            "[5, 3]": 1,
            "[4, 1, 3]": 1,
            "[4, 1, 1, 2, 5]": 1,
            "[8, 1, 2, 5]": 1,
            "[8, 3]": 1,
            "[6, 2, 3, 2]": 1,
            "[3, 3, 5, 2]": 1,
            "[3, 3, 1, 4, 2]": 1
        })
        self.assertEqual(len(results.psc_duration), 9)
        self.assertAlmostEqual(results.psc_duration["[2]"], Decimal(0.5))
        self.assertAlmostEqual(results.psc_duration["[5, 3]"], Decimal(0.625))
        self.assertAlmostEqual(results.psc_duration["[4, 1, 3]"], Decimal(1/3-0.125))
        self.assertAlmostEqual(results.psc_duration["[4, 1, 1, 2, 5]"], Decimal(1/6))
        self.assertAlmostEqual(results.psc_duration["[8, 1, 2, 5]"], Decimal(1/6))
        self.assertAlmostEqual(results.psc_duration["[8, 3]"], Decimal(1/3+0.4))
        self.assertAlmostEqual(results.psc_duration["[6, 2, 3, 2]"], Decimal(0.1))
        self.assertAlmostEqual(results.psc_duration["[3, 3, 5, 2]"], Decimal(0.3))
        self.assertAlmostEqual(results.psc_duration["[3, 3, 1, 4, 2]"], Decimal(6/5))
        self.assertEqual(results.chord_spacing_contour_frequency, {
            "<0>": 1,
            "<1, 0>": 2,
            "<2, 0, 1>": 1,
            "<2, 0, 0, 1, 3>": 1,
            "<3, 0, 1, 2>": 1,
            "<2, 0, 1, 0>": 1,
            "<1, 1, 2, 0>": 1,
            "<2, 2, 0, 3, 1>": 1
        })
        self.assertEqual(len(results.chord_spacing_contour_duration), 8)
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<0>"], Decimal(0.5))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<1, 0>"], Decimal(1.025+1/3))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<2, 0, 1>"], Decimal(1/3-0.125))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<2, 0, 0, 1, 3>"], Decimal(1/6))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<3, 0, 1, 2>"], Decimal(1/6))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<2, 0, 1, 0>"], Decimal(0.1))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<1, 1, 2, 0>"], Decimal(0.3))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<2, 2, 0, 3, 1>"], Decimal(6/5))

    def test_carter2(self):
        """
        Tests measure 130 of Carter V
        """
        results = salami_slice_analyze.analyze(Path(__file__).parent / "data/test_carter2.musicxml")[0]
        slices = results.slices
        self.assertEqual(len(slices), 9)

        self.assertEqual(slices[0].pseg, pseg.make_pseg12(3, 7))
        self.assertEqual(slices[0].pset, pset.make_pset12(3, 7))
        self.assertEqual(slices[0].psets, 
            [set(), set(), pset.make_pset12(3, 7), set()])
        self.assertEqual(slices[0].pcset, pcset.make_pcset12(3, 7))
        self.assertEqual(slices[0].pcseg, pcseg.make_pcseg12(3, 7))
        self.assertEqual(slices[0].pitchseg, [3, 7])
        self.assertEqual(slices[0].sc_name, "(2-4)[04]")
        self.assertEqual(slices[0].chord_spacing_contour, [0])
        self.assertTrue(np.isnan(slices[0].chord_spacing_index))

        self.assertEqual(slices[1].pseg, pseg.make_pseg12(3, 7, 13))
        self.assertEqual(slices[1].pset, pset.make_pset12(3, 7, 13))
        self.assertEqual(slices[1].psets, 
            [set(), {Pitch(13),}, pset.make_pset12(3, 7), set()])
        self.assertEqual(slices[1].pcset, pcset.make_pcset12(1, 3, 7))
        self.assertEqual(slices[1].pcseg, pcseg.make_pcseg12(3, 7, 1))
        self.assertEqual(slices[1].pitchseg, [3, 7, 13])
        self.assertEqual(slices[1].sc_name, "(3-8)[026]")
        self.assertEqual(slices[1].chord_spacing_contour, [0, 1])
        self.assertAlmostEqual(slices[1].chord_spacing_index, csi(slices[1].pseg))
        
        self.assertEqual(slices[2].pseg, pseg.make_pseg12(13))
        self.assertEqual(slices[2].pset, pset.make_pset12(13))
        self.assertEqual(slices[2].psets, [set(), {Pitch(13),}, set(), set()])
        self.assertEqual(slices[2].pcset, {PitchClass(1),})
        self.assertEqual(slices[2].pcseg, [PitchClass(1),])
        self.assertEqual(slices[2].pitchseg, [13,])
        self.assertEqual(slices[2].sc_name, "(1-1)[0]")
        self.assertEqual(slices[2].chord_spacing_contour, [])
        self.assertTrue(np.isnan(slices[2].chord_spacing_index))
        
        self.assertEqual(slices[3].pseg, pseg.make_pseg12(12))
        self.assertEqual(slices[3].pset, pset.make_pset12(12))
        self.assertEqual(slices[3].psets, [set(), {Pitch(12),}, set(), set()])
        self.assertEqual(slices[3].pcset, {PitchClass(0),})
        self.assertEqual(slices[3].pcseg, [PitchClass(0),])
        self.assertEqual(slices[3].pitchseg, [12,])
        self.assertEqual(slices[3].sc_name, "(1-1)[0]")
        self.assertEqual(slices[3].chord_spacing_contour, [])
        self.assertTrue(np.isnan(slices[3].chord_spacing_index))
        
        self.assertEqual(slices[4].pseg, [])
        self.assertEqual(slices[4].pset, set())
        self.assertEqual(slices[4].psets, [set(), set(), set(), set()])
        self.assertEqual(slices[4].pcset, set())
        self.assertEqual(slices[4].pcseg, [])
        self.assertEqual(slices[4].pitchseg, [])
        self.assertEqual(slices[4].sc_name, "(0-1)[]")
        self.assertEqual(slices[4].chord_spacing_contour, [])
        self.assertTrue(np.isnan(slices[4].chord_spacing_index))
        
        self.assertEqual(slices[5].pseg, pseg.make_pseg12(-4, 5, 14))
        self.assertEqual(slices[5].pset, pset.make_pset12(-4, 5, 14))
        self.assertEqual(slices[5].psets, 
            [set(), set(), pset.make_pset12(-4, 5, 14), set()])
        self.assertEqual(slices[5].pcset, pcset.make_pcset12(2, 5, 8))
        self.assertEqual(slices[5].pcseg, pcseg.make_pcseg12(8, 5, 2))
        self.assertEqual(slices[5].pitchseg, [-4, 5, 14])
        self.assertEqual(slices[5].sc_name, "(3-10)[036]")
        self.assertEqual(slices[5].chord_spacing_contour, [0, 0])
        self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))
        
        self.assertEqual(slices[6].pseg, [])
        self.assertEqual(slices[6].pset, set())
        self.assertEqual(slices[6].psets, [set(), set(), set(), set()])
        self.assertEqual(slices[6].pcset, set())
        self.assertEqual(slices[6].pcseg, [])
        self.assertEqual(slices[6].pitchseg, [])
        self.assertEqual(slices[6].sc_name, "(0-1)[]")
        self.assertEqual(slices[6].chord_spacing_contour, [])
        self.assertTrue(np.isnan(slices[6].chord_spacing_index))
        
        self.assertEqual(slices[7].pseg, pseg.make_pseg12(-5, 12, 15, 22))
        self.assertEqual(slices[7].pset, pset.make_pset12(-5, 12, 15, 22))
        self.assertEqual(slices[7].psets, [set(), pset.make_pset12(-5, 12, 15, 22), set(), set()])
        self.assertEqual(slices[7].pcset, pcset.make_pcset12(0, 3, 7, 10))
        self.assertEqual(slices[7].pcseg, pcseg.make_pcseg12(7, 0, 3, 10))
        self.assertEqual(slices[7].pitchseg, [-5, 12, 15, 22])
        self.assertEqual(slices[7].sc_name, "(4-26)[0358]")
        self.assertEqual(slices[7].chord_spacing_contour, [2, 0, 1])
        self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))
        
        self.assertEqual(slices[8].pseg, [Pitch(22),])
        self.assertEqual(slices[8].pset, {Pitch(22),})
        self.assertEqual(slices[8].psets, [set(), {Pitch(22),}, set(), set()])
        self.assertEqual(slices[8].pcset, {PitchClass(10),})
        self.assertEqual(slices[8].pcseg, [PitchClass(10),])
        self.assertEqual(slices[8].pitchseg, [22,])
        self.assertEqual(slices[8].sc_name, "(1-1)[0]")
        self.assertEqual(slices[8].chord_spacing_contour, [])
        self.assertTrue(np.isnan(slices[8].chord_spacing_index))

        self.assertAlmostEqual(slices[0].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[1].duration, Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(slices[2].duration, Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(slices[3].duration, Decimal(60/115.2))
        self.assertAlmostEqual(slices[4].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[5].duration, Decimal(60/115.2))
        self.assertAlmostEqual(slices[6].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[7].duration, Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(slices[8].duration, Decimal((1/3) * (60/115.2)))

        self.assertEqual(slices[0].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[1].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[2].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[3].quarter_duration, Fraction(1, 1))
        self.assertEqual(slices[4].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[5].quarter_duration, Fraction(1, 1))
        self.assertEqual(slices[6].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[7].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[8].quarter_duration, Fraction(1, 3))

        self.assertAlmostEqual(results.chord_spacing_index_avg, (
            float(slices[1].duration) * csi(slices[1].pseg) + 
            float(slices[5].duration) * csi(slices[5].pseg) + 
            float(slices[7].duration) * csi(slices[7].pseg)
        ) / (float(slices[1].duration) + float(slices[5].duration) + float(slices[7].duration)))
        self.assertAlmostEqual(float(results.duration_avg), 
            np.average([float(s.duration) for s in slices]))
        self.assertAlmostEqual(float(results.duration), 60/115.2*4)
        self.assertEqual(results.quarter_duration, 4)
        self.assertEqual(results.measure_num_first, 1)
        self.assertEqual(results.measure_num_last, 1)
        self.assertEqual(results.num_measures, 1)
        self.assertEqual(results.num_voices, 4)
        self.assertEqual(results.max_pitch_count_with_duplicates, 4)
        self.assertEqual(results.lower_bound, -5)
        self.assertEqual(results.upper_bound, 22)
        self.assertEqual(results.start_time, 0)
        self.assertEqual(len(results.pc_duration), 8)
        self.assertAlmostEqual(results.pc_duration[0], Decimal((4/3) * (60/115.2)))
        self.assertAlmostEqual(results.pc_duration[1], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.pc_duration[2], Decimal(60/115.2))
        self.assertAlmostEqual(results.pc_duration[3], Decimal((5/6) * (60/115.2)))
        self.assertAlmostEqual(results.pc_duration[5], Decimal(60/115.2))
        self.assertAlmostEqual(results.pc_duration[7], Decimal((5/6) * (60/115.2)))
        self.assertAlmostEqual(results.pc_duration[8], Decimal(60/115.2))
        self.assertAlmostEqual(results.pc_duration[10], Decimal((2/3) * (60/115.2)))
        self.assertEqual(len(results.pitch_duration), 10)
        self.assertAlmostEqual(results.pitch_duration[-5], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[-4], Decimal(60/115.2))
        self.assertAlmostEqual(results.pitch_duration[3], Decimal(0.5 * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[5], Decimal(60/115.2))
        self.assertAlmostEqual(results.pitch_duration[7], Decimal(0.5 * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[12], Decimal((4/3) * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[13], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[14], Decimal(60/115.2))
        self.assertAlmostEqual(results.pitch_duration[15], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.pitch_duration[22], Decimal((2/3) * (60/115.2)))
        self.assertEqual(results.pc_frequency, {
            0: 2, 1: 1, 2: 1, 3: 2, 5: 1, 7: 2, 8: 1, 10: 1
        })
        self.assertEqual(results.pitch_frequency, {
            -5: 1, -4: 1, 3: 1, 5: 1, 7: 1, 12: 2, 13: 1, 14: 1, 15: 1, 22: 1
        })
        self.assertEqual(results.pcsc_frequency, {
            "(2-4)[04]": 1,
            "(3-8)[026]": 1,
            "(1-1)[0]": 3,
            "(0-1)[]": 2,
            "(3-10)[036]": 1,
            "(4-26)[0358]": 1
        })
        self.assertEqual(len(results.pcsc_duration), 6)
        self.assertAlmostEqual(results.pcsc_duration["(2-4)[04]"], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.pcsc_duration["(3-8)[026]"], Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(results.pcsc_duration["(1-1)[0]"], Decimal(1.5 * (60/115.2)))
        self.assertAlmostEqual(results.pcsc_duration["(0-1)[]"], Decimal((2/3) * (60/115.2)))
        self.assertAlmostEqual(results.pcsc_duration["(3-10)[036]"], Decimal(60/115.2))
        self.assertAlmostEqual(results.pcsc_duration["(4-26)[0358]"], Decimal((1/3) * (60/115.2)))
        self.assertEqual(results.psc_frequency, {
            "[4]": 1,
            "[4, 6]": 1,
            "[]": 5,
            "[9, 9]": 1,
            "[17, 3, 7]": 1
        })
        self.assertEqual(len(results.psc_duration), 5)
        self.assertAlmostEqual(results.psc_duration["[4]"], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.psc_duration["[4, 6]"], Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(results.psc_duration["[]"], Decimal((13/6) * (60/115.2)))
        self.assertAlmostEqual(results.psc_duration["[9, 9]"], Decimal(60/115.2))
        self.assertAlmostEqual(results.psc_duration["[17, 3, 7]"], Decimal((1/3) * (60/115.2)))
        self.assertEqual(results.chord_spacing_contour_frequency, {
            "<0>": 1,
            "<0, 1>": 1,
            "<>": 5,
            "<0, 0>": 1,
            "<2, 0, 1>": 1
        })
        self.assertEqual(len(results.chord_spacing_contour_duration), 5)
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<0>"], Decimal((1/3) * (60/115.2)))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<0, 1>"], Decimal((1/6) * (60/115.2)))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<>"], Decimal((13/6) * (60/115.2)))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<0, 0>"], Decimal(60/115.2))
        self.assertAlmostEqual(results.chord_spacing_contour_duration["<2, 0, 1>"], Decimal((1/3) * (60/115.2)))
        
    def test_carter3(self):
        """
        Tests measures 101-103 of Carter V
        """
        results = salami_slice_analyze.analyze(Path(__file__).parent / "data/test_carter3.musicxml")[0]
        slices = results.slices
        self.assertEqual(len(slices), 18)

        self.assertEqual(slices[0].pseg, pseg.make_pseg12(-20, -13, -7, 0, 2, 10, 18, 32))
        self.assertEqual(slices[0].pset, pset.make_pset12(-20, -13, -7, 0, 2, 10, 18, 32))
        self.assertEqual(slices[0].psets, 
            [pset.make_pset12(18, 32), pset.make_pset12(0, 10), 
            pset.make_pset12(-7, 2), pset.make_pset12(-20, -13)])
        self.assertEqual(slices[0].pcset, pcset.make_pcset12(0, 2, 4, 5, 6, 8, 10, 11))
        self.assertEqual(slices[0].pcseg, pcseg.make_pcseg12(4, 11, 5, 0, 2, 10, 6, 8))
        self.assertEqual(slices[0].pitchseg, [-20, -13, -7, 0, 2, 10, 18, 32])
        self.assertEqual(slices[0].sc_name, "(8-25)[0124678A]")
        self.assertEqual(slices[0].chord_spacing_contour, [2, 1, 2, 0, 3, 3, 4])
        # psc 7 6 7 2 8 8 14
        self.assertAlmostEqual(slices[0].chord_spacing_index, csi(slices[0].pseg))

        self.assertEqual(slices[1].pseg, [])
        self.assertEqual(slices[1].pset, set())
        self.assertEqual(slices[1].psets, [set(), set(), set(), set()])
        self.assertEqual(slices[1].pcset, set())
        self.assertEqual(slices[1].pcseg, [])
        self.assertEqual(slices[1].pitchseg, [])
        self.assertEqual(slices[1].sc_name, "(0-1)[]")
        self.assertEqual(slices[1].chord_spacing_contour, [])
        # psc []
        self.assertTrue(np.isnan(slices[1].chord_spacing_index))

        self.assertEqual(slices[2].pseg, [Pitch(32)])
        self.assertEqual(slices[2].pset, {Pitch(32),})
        self.assertEqual(slices[2].psets, 
            [{Pitch(32),}, set(), set(), set()])
        self.assertEqual(slices[2].pcset, {PitchClass(8),})
        self.assertEqual(slices[2].pcseg, [PitchClass(8)])
        self.assertEqual(slices[2].pitchseg, [32])
        self.assertEqual(slices[2].sc_name, "(1-1)[0]")
        self.assertEqual(slices[2].chord_spacing_contour, [])
        # psc []
        self.assertTrue(np.isnan(slices[2].chord_spacing_index))

        self.assertEqual(slices[3].pseg, pseg.make_pseg12(-7, 32))
        self.assertEqual(slices[3].pset, pset.make_pset12(-7, 32))
        self.assertEqual(slices[3].psets, 
            [{Pitch(32),}, set(), {Pitch(-7),}, set()])
        self.assertEqual(slices[3].pcset, pcset.make_pcset12(5, 8))
        self.assertEqual(slices[3].pcseg, pcseg.make_pcseg12(5, 8))
        self.assertEqual(slices[3].pitchseg, [-7, 32])
        self.assertEqual(slices[3].sc_name, "(2-3)[03]")
        self.assertEqual(slices[3].chord_spacing_contour, [0])
        # psc [39]
        self.assertTrue(np.isnan(slices[3].chord_spacing_index))

        self.assertEqual(slices[4].pseg, pseg.make_pseg12(-7, 0, 32))
        self.assertEqual(slices[4].pset, pset.make_pset12(-7, 0, 32))
        self.assertEqual(slices[4].psets, 
            [{Pitch(32),}, {Pitch(0),}, {Pitch(-7),}, set()])
        self.assertEqual(slices[4].pcset, pcset.make_pcset12(0, 5, 8))
        self.assertEqual(slices[4].pcseg, pcseg.make_pcseg12(5, 0, 8))
        self.assertEqual(slices[4].pitchseg, [-7, 0, 32])
        self.assertEqual(slices[4].sc_name, "(3-11)[037]")
        self.assertEqual(slices[4].chord_spacing_contour, [0, 1])
        # psc [7, 32]
        self.assertAlmostEqual(slices[4].chord_spacing_index, csi(slices[4].pseg))

        self.assertEqual(slices[5].pseg, pseg.make_pseg12(-20, -13, -7, 0, 32))
        self.assertEqual(slices[5].pset, pset.make_pset12(-20, -13, -7, 0, 32))
        self.assertEqual(slices[5].psets, 
            [{Pitch(32),}, {Pitch(0),}, {Pitch(-7),}, pset.make_pset12(-20, -13)])
        self.assertEqual(slices[5].pcset, pcset.make_pcset12(0, 4, 5, 8, 11))
        self.assertEqual(slices[5].pcseg, pcseg.make_pcseg12(4, 11, 5, 0, 8))
        self.assertEqual(slices[5].pitchseg, [-20, -13, -7, 0, 32])
        self.assertEqual(slices[5].sc_name, "(5-22)[01478]")
        self.assertEqual(slices[5].chord_spacing_contour, [1, 0, 1, 2])
        # psc [7, 6, 7, 32]
        self.assertAlmostEqual(slices[5].chord_spacing_index, csi(slices[5].pseg))

        self.assertEqual(slices[6].pseg, pseg.make_pseg12(-20, -13, -7, 32))
        self.assertEqual(slices[6].pset, pset.make_pset12(-20, -13, -7, 32))
        self.assertEqual(slices[6].psets, 
            [{Pitch(32),}, set(), {Pitch(-7),}, pset.make_pset12(-20, -13)])
        self.assertEqual(slices[6].pcset, pcset.make_pcset12(4, 5, 8, 11))
        self.assertEqual(slices[6].pcseg, pcseg.make_pcseg12(4, 11, 5, 8))
        self.assertEqual(slices[6].pitchseg, [-20, -13, -7, 32])
        self.assertEqual(slices[6].sc_name, "(4-18)[0147]")
        self.assertEqual(slices[6].chord_spacing_contour, [1, 0, 2])
        # psc [7, 6, 39]
        self.assertAlmostEqual(slices[6].chord_spacing_index, csi(slices[6].pseg))

        self.assertEqual(slices[7].pseg, pseg.make_pseg12(-20, -13, 32))
        self.assertEqual(slices[7].pset, pset.make_pset12(-20, -13, 32))
        self.assertEqual(slices[7].psets, 
            [{Pitch(32),}, set(), set(), pset.make_pset12(-20, -13)])
        self.assertEqual(slices[7].pcset, pcset.make_pcset12(4, 8, 11))
        self.assertEqual(slices[7].pcseg, pcseg.make_pcseg12(4, 11, 8))
        self.assertEqual(slices[7].pitchseg, [-20, -13, 32])
        self.assertEqual(slices[7].sc_name, "(3-11)[037]")
        self.assertEqual(slices[7].chord_spacing_contour, [0, 1])
        # psc [7, 45]
        self.assertAlmostEqual(slices[7].chord_spacing_index, csi(slices[7].pseg))

        self.assertEqual(slices[8].pseg, pseg.make_pseg12(-20, -13, 18))
        self.assertEqual(slices[8].pset, pset.make_pset12(-20, -13, 18))
        self.assertEqual(slices[8].psets, 
            [{Pitch(18),}, set(), set(), pset.make_pset12(-20, -13)])
        self.assertEqual(slices[8].pcset, pcset.make_pcset12(4, 6, 11))
        self.assertEqual(slices[8].pcseg, pcseg.make_pcseg12(4, 11, 6))
        self.assertEqual(slices[8].pitchseg, [-20, -13, 18])
        self.assertEqual(slices[8].sc_name, "(3-9)[027]")
        self.assertEqual(slices[8].chord_spacing_contour, [0, 1])
        # psc [7, 45]
        self.assertAlmostEqual(slices[8].chord_spacing_index, csi(slices[8].pseg))

        self.assertEqual(slices[9].pseg, pseg.make_pseg12(-20, -13, -3, 10, 18))
        self.assertEqual(slices[9].pset, pset.make_pset12(-20, -13, -3, 10, 18))
        self.assertEqual(slices[9].psets, 
            [{Pitch(18),}, pset.make_pset12(-3, 10), set(), pset.make_pset12(-20, -13)])
        self.assertEqual(slices[9].pcset, pcset.make_pcset12(4, 6, 9, 10, 11))
        self.assertEqual(slices[9].pcseg, pcseg.make_pcseg12(4, 11, 9, 10, 6))
        self.assertEqual(slices[9].pitchseg, [-20, -13, -3, 10, 18])
        self.assertEqual(slices[9].sc_name, "(5-14)[01257]")
        self.assertEqual(slices[9].chord_spacing_contour, [0, 2, 3, 1])
        # psc [7, 10, 13, 8]
        self.assertAlmostEqual(slices[9].chord_spacing_index, csi(slices[9].pseg))

        self.assertEqual(slices[10].pseg, pseg.make_pseg12(-20, -13, -7, -3, 10, 18))
        self.assertEqual(slices[10].pset, pset.make_pset12(-20, -13, -7, -3, 10, 18))
        self.assertEqual(slices[10].psets, 
            [{Pitch(18),}, pset.make_pset12(-3, 10), {Pitch(-7),}, pset.make_pset12(-20, -13)])
        self.assertEqual(slices[10].pcset, pcset.make_pcset12(4, 5, 6, 9, 10, 11))
        self.assertEqual(slices[10].pcseg, pcseg.make_pcseg12(4, 11, 5, 9, 10, 6))
        self.assertEqual(slices[10].pitchseg, [-20, -13, -7, -3, 10, 18])
        self.assertEqual(slices[10].sc_name, "(6-Z6)[012567]")
        self.assertEqual(slices[10].chord_spacing_contour, [2, 1, 0, 4, 3])
        # psc [7, 6, 4, 13, 8]
        self.assertAlmostEqual(slices[10].chord_spacing_index, csi(slices[10].pseg))

        self.assertEqual(slices[11].pseg, pseg.make_pseg12(-20, -13, -7, -3, 10))
        self.assertEqual(slices[11].pset, pset.make_pset12(-20, -13, -7, -3, 10))
        self.assertEqual(slices[11].psets, 
            [set(), pset.make_pset12(-3, 10), {Pitch(-7),}, pset.make_pset12(-20, -13)])
        self.assertEqual(slices[11].pcset, pcset.make_pcset12(4, 5, 9, 10, 11))
        self.assertEqual(slices[11].pcseg, pcseg.make_pcseg12(4, 11, 5, 9, 10))
        self.assertEqual(slices[11].pitchseg, [-20, -13, -7, -3, 10])
        self.assertEqual(slices[11].sc_name, "(5-7)[01267]")
        self.assertEqual(slices[11].chord_spacing_contour, [2, 1, 0, 3])
        # psc [7, 6, 4, 13]
        self.assertAlmostEqual(slices[11].chord_spacing_index, csi(slices[11].pseg))

        self.assertEqual(slices[12].pseg, pseg.make_pseg12(-7, -3, 10))
        self.assertEqual(slices[12].pset, pset.make_pset12(-7, -3, 10))
        self.assertEqual(slices[12].psets, 
            [set(), pset.make_pset12(-3, 10), {Pitch(-7),}, set()])
        self.assertEqual(slices[12].pcset, pcset.make_pcset12(5, 9, 10))
        self.assertEqual(slices[12].pcseg, pcseg.make_pcseg12(5, 9, 10))
        self.assertEqual(slices[12].pitchseg, [-7, -3, 10])
        self.assertEqual(slices[12].sc_name, "(3-4)[015]")
        self.assertEqual(slices[12].chord_spacing_contour, [0, 1])
        # psc [4, 13]
        self.assertAlmostEqual(slices[12].chord_spacing_index, csi(slices[12].pseg))

        self.assertEqual(slices[13].pseg, pseg.make_pseg12(-7, -3, 10))
        self.assertEqual(slices[13].pset, pset.make_pset12(-7, -3, 10))
        self.assertEqual(slices[13].psets, 
            [set(), pset.make_pset12(-3, 10), pset.make_pset12(-7, -3), set()])
        self.assertEqual(slices[13].pcset, pcset.make_pcset12(5, 9, 10))
        self.assertEqual(slices[13].pcseg, pcseg.make_pcseg12(5, 9, 10))
        self.assertEqual(slices[13].pitchseg, [-7, -3, -3, 10])
        self.assertEqual(slices[13].sc_name, "(3-4)[015]")
        self.assertEqual(slices[13].chord_spacing_contour, [0, 1])
        # psc [4, 13]
        self.assertAlmostEqual(slices[13].chord_spacing_index, csi(slices[13].pseg))

        self.assertEqual(slices[14].pseg, pseg.make_pseg12(-7, -3))
        self.assertEqual(slices[14].pset, pset.make_pset12(-7, -3))
        self.assertEqual(slices[14].psets, 
            [set(), set(), pset.make_pset12(-7, -3), set()])
        self.assertEqual(slices[14].pcset, pcset.make_pcset12(5, 9))
        self.assertEqual(slices[14].pcseg, pcseg.make_pcseg12(5, 9))
        self.assertEqual(slices[14].pitchseg, [-7, -3])
        self.assertEqual(slices[14].sc_name, "(2-4)[04]")
        self.assertEqual(slices[14].chord_spacing_contour, [0])
        # psc [4]
        self.assertTrue(np.isnan(slices[14].chord_spacing_index))

        self.assertEqual(slices[15].pseg, pseg.make_pseg12(-7, -3, 0, 7))
        self.assertEqual(slices[15].pset, pset.make_pset12(-7, -3, 0, 7))
        self.assertEqual(slices[15].psets, 
            [set(), pset.make_pset12(0, 7), pset.make_pset12(-7, -3), set()])
        self.assertEqual(slices[15].pcset, pcset.make_pcset12(0, 5, 7, 9))
        self.assertEqual(slices[15].pcseg, pcseg.make_pcseg12(5, 9, 0, 7))
        self.assertEqual(slices[15].pitchseg, [-7, -3, 0, 7])
        self.assertEqual(slices[15].sc_name, "(4-22)[0247]")
        self.assertEqual(slices[15].chord_spacing_contour, [1, 0, 2])
        # psc [4, 3, 7]
        self.assertAlmostEqual(slices[15].chord_spacing_index, csi(slices[15].pseg))

        self.assertEqual(slices[16].pseg, pseg.make_pseg12(-13, -11, -7, -3, 0, 7))
        self.assertEqual(slices[16].pset, pset.make_pset12(-13, -11, -7, -3, 0, 7))
        self.assertEqual(slices[16].psets, 
            [set(), pset.make_pset12(0, 7), pset.make_pset12(-7, -3), pset.make_pset12(-13, -11)])
        self.assertEqual(slices[16].pcset, pcset.make_pcset12(0, 1, 5, 7, 9, 11))
        self.assertEqual(slices[16].pcseg, pcseg.make_pcseg12(11, 1, 5, 9, 0, 7))
        self.assertEqual(slices[16].pitchseg, [-13, -11, -7, -3, 0, 7])
        self.assertEqual(slices[16].sc_name, "(6-22)[012468]")
        self.assertEqual(slices[16].chord_spacing_contour, [0, 2, 2, 1, 3])
        # psc [2, 4, 4, 3, 7]
        self.assertAlmostEqual(slices[16].chord_spacing_index, csi(slices[16].pseg))

        self.assertEqual(slices[17].pseg, pseg.make_pseg12(-13, -11, -7, -3, 0, 7, 15, 18))
        self.assertEqual(slices[17].pset, pset.make_pset12(-13, -11, -7, -3, 0, 7, 15, 18))
        self.assertEqual(slices[17].psets, 
            [pset.make_pset12(15, 18), pset.make_pset12(0, 7), pset.make_pset12(-7, -3), pset.make_pset12(-13, -11)])
        self.assertEqual(slices[17].pcset, pcset.make_pcset12(0, 1, 3, 5, 6, 7, 9, 11))
        self.assertEqual(slices[17].pcseg, pcseg.make_pcseg12(11, 1, 5, 9, 0, 7, 3, 6))
        self.assertEqual(slices[17].pitchseg, [-13, -11, -7, -3, 0, 7, 15, 18])
        self.assertEqual(slices[17].sc_name, "(8-25)[0124678A]")
        self.assertEqual(slices[17].chord_spacing_contour, [0, 2, 2, 1, 3, 4, 1])
        # psc [2, 4, 4, 3, 7, 8, 3]
        self.assertAlmostEqual(slices[17].chord_spacing_index, csi(slices[17].pseg))

        self.assertAlmostEqual(slices[0].duration, Decimal(2))
        self.assertAlmostEqual(slices[1].duration, Decimal(1/3))
        self.assertAlmostEqual(slices[2].duration, Decimal(11/21))
        self.assertAlmostEqual(slices[3].duration, Decimal(19/35))
        self.assertAlmostEqual(slices[4].duration, Decimal(17/20))
        self.assertAlmostEqual(slices[5].duration, Decimal(3/20))
        self.assertAlmostEqual(slices[6].duration, Decimal(1/35))
        self.assertAlmostEqual(slices[7].duration, Decimal(5/21))
        self.assertAlmostEqual(slices[8].duration, Decimal(4/3))
        self.assertAlmostEqual(slices[9].duration, Decimal(1/7))
        self.assertAlmostEqual(slices[10].duration, Decimal(4/21))
        self.assertAlmostEqual(slices[11].duration, Decimal(1/6))
        self.assertAlmostEqual(slices[12].duration, Decimal(15/14))
        self.assertAlmostEqual(slices[13].duration, Decimal(8/35))
        self.assertAlmostEqual(slices[14].duration, Decimal(1))
        self.assertAlmostEqual(slices[15].duration, Decimal(6/5))
        self.assertAlmostEqual(slices[16].duration, Decimal(4/3))
        self.assertAlmostEqual(slices[17].duration, Decimal(2/3))

        self.assertEqual(slices[0].quarter_duration, Fraction(2, 1))
        self.assertEqual(slices[1].quarter_duration, Fraction(1, 3))
        self.assertEqual(slices[2].quarter_duration, Fraction(11, 21))
        self.assertEqual(slices[3].quarter_duration, Fraction(19, 35))
        self.assertEqual(slices[4].quarter_duration, Fraction(17, 20))
        self.assertEqual(slices[5].quarter_duration, Fraction(3, 20))
        self.assertEqual(slices[6].quarter_duration, Fraction(1, 35))
        self.assertEqual(slices[7].quarter_duration, Fraction(5, 21))
        self.assertEqual(slices[8].quarter_duration, Fraction(4, 3))
        self.assertEqual(slices[9].quarter_duration, Fraction(1, 7))
        self.assertEqual(slices[10].quarter_duration, Fraction(4, 21))
        self.assertEqual(slices[11].quarter_duration, Fraction(1, 6))
        self.assertEqual(slices[12].quarter_duration, Fraction(15, 14))
        self.assertEqual(slices[13].quarter_duration, Fraction(8, 35))
        self.assertEqual(slices[14].quarter_duration, Fraction(1, 1))
        self.assertEqual(slices[15].quarter_duration, Fraction(6, 5))
        self.assertEqual(slices[16].quarter_duration, Fraction(4, 3))
        self.assertEqual(slices[17].quarter_duration, Fraction(2, 3))

        non_nan_csi = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17]
        self.assertAlmostEqual(results.chord_spacing_index_avg, 
            sum([float(slices[i].duration) * csi(slices[i].pseg) for i in non_nan_csi]) / 
            sum([float(slices[i].duration) for i in non_nan_csi]))
        self.assertAlmostEqual(float(results.duration_avg), 
            np.average([float(s.duration) for s in slices]))
        self.assertAlmostEqual(float(results.duration), 12.0)
        self.assertEqual(results.quarter_duration, 12)
        self.assertEqual(results.measure_num_first, 1)
        self.assertEqual(results.measure_num_last, 3)
        self.assertEqual(results.num_measures, 3)
        self.assertEqual(results.num_voices, 4)
        self.assertEqual(results.max_pitch_count_with_duplicates, 8)
        self.assertEqual(results.lower_bound, -20)
        self.assertEqual(results.upper_bound, 32)
        self.assertEqual(results.start_time, 0)
        self.assertEqual(len(results.pc_duration), 12)
        self.assertAlmostEqual(results.pc_duration[0], Decimal(31/5))
        self.assertAlmostEqual(results.pc_duration[1], Decimal(2))
        self.assertAlmostEqual(results.pc_duration[2], Decimal(2))
        self.assertAlmostEqual(results.pc_duration[3], Decimal(2/3))
        self.assertAlmostEqual(results.pc_duration[4], Decimal(17/4))
        self.assertAlmostEqual(results.pc_duration[5], Decimal(66/7))
        self.assertAlmostEqual(results.pc_duration[6], Decimal(13/3))
        self.assertAlmostEqual(results.pc_duration[7], Decimal(16/5))
        self.assertAlmostEqual(results.pc_duration[8], Decimal(13/3))
        self.assertAlmostEqual(results.pc_duration[9], Decimal(6))
        self.assertAlmostEqual(results.pc_duration[10], Decimal(19/5))
        self.assertAlmostEqual(results.pc_duration[11], Decimal(25/4))
        self.assertEqual(len(results.pitch_duration), 12)
        self.assertAlmostEqual(results.pitch_duration[0], Decimal(31/5))
        self.assertAlmostEqual(results.pitch_duration[-11], Decimal(2))
        self.assertAlmostEqual(results.pitch_duration[2], Decimal(2))
        self.assertAlmostEqual(results.pitch_duration[15], Decimal(2/3))
        self.assertAlmostEqual(results.pitch_duration[-20], Decimal(17/4))
        self.assertAlmostEqual(results.pitch_duration[-7], Decimal(66/7))
        self.assertAlmostEqual(results.pitch_duration[18], Decimal(13/3))
        self.assertAlmostEqual(results.pitch_duration[7], Decimal(16/5))
        self.assertAlmostEqual(results.pitch_duration[32], Decimal(13/3))
        self.assertAlmostEqual(results.pitch_duration[-3], Decimal(6))
        self.assertAlmostEqual(results.pitch_duration[10], Decimal(19/5))
        self.assertAlmostEqual(results.pitch_duration[-13], Decimal(25/4))
        self.assertEqual(results.pc_frequency, {
            0: 3, 1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 3, 7: 1, 8: 2, 9: 1, 10: 2, 11: 3
        })
        self.assertEqual(results.pitch_frequency, {
            -20: 2, -13: 3, -11: 1, -7: 3, -3: 1, 0: 3, 2: 1, 7: 1, 10: 2, 15: 1, 18: 3, 32: 2
        })
        self.assertEqual(results.pcsc_frequency, {
            "(8-25)[0124678A]": 2, #0, 17
            "(0-1)[]": 1, #1
            "(1-1)[0]": 1, #2
            "(2-3)[03]": 1, #3
            "(3-11)[037]": 2, #4, 7
            "(5-22)[01478]": 1, #5
            "(4-18)[0147]": 1, #6
            "(3-9)[027]": 1, #8
            "(5-14)[01257]": 1, #9
            "(6-Z6)[012567]": 1, #10
            "(5-7)[01267]": 1, #11
            "(3-4)[015]": 2, #12, 13
            "(2-4)[04]": 1, #14
            "(4-22)[0247]": 1, #15
            "(6-22)[012468]": 1, #16
        })
        self.assertEqual(len(results.pcsc_duration), 15)
        self.assertAlmostEqual(results.pcsc_duration["(8-25)[0124678A]"], Decimal(8/3))
        self.assertAlmostEqual(results.pcsc_duration["(0-1)[]"], Decimal(1/3))
        self.assertAlmostEqual(results.pcsc_duration["(1-1)[0]"], Decimal(11/21))
        self.assertAlmostEqual(results.pcsc_duration["(2-3)[03]"], Decimal(19/35))
        self.assertAlmostEqual(results.pcsc_duration["(3-11)[037]"], Decimal(457/420))
        self.assertAlmostEqual(results.pcsc_duration["(5-22)[01478]"], Decimal(3/20))
        self.assertAlmostEqual(results.pcsc_duration["(4-18)[0147]"], Decimal(1/35))
        self.assertAlmostEqual(results.pcsc_duration["(3-9)[027]"], Decimal(4/3))
        self.assertAlmostEqual(results.pcsc_duration["(5-14)[01257]"], Decimal(1/7))
        self.assertAlmostEqual(results.pcsc_duration["(6-Z6)[012567]"], Decimal(4/21))
        self.assertAlmostEqual(results.pcsc_duration["(5-7)[01267]"], Decimal(1/6))
        self.assertAlmostEqual(results.pcsc_duration["(3-4)[015]"], Decimal(13/10))
        self.assertAlmostEqual(results.pcsc_duration["(2-4)[04]"], Decimal(1))
        self.assertAlmostEqual(results.pcsc_duration["(4-22)[0247]"], Decimal(6/5))
        self.assertAlmostEqual(results.pcsc_duration["(6-22)[012468]"], Decimal(4/3))
        # self.assertEqual(results.psc_frequency, {
        #     "[4]": 1,
        #     "[4, 6]": 1,
        #     "[]": 5,
        #     "[9, 9]": 1,
        #     "[17, 3, 7]": 1
        # })
        # self.assertEqual(len(results.psc_duration), 5)
        # self.assertAlmostEqual(results.psc_duration["[4]"], Decimal((1/3) * (60/115.2)))
        
if __name__ == "__main__":
    unittest.main()