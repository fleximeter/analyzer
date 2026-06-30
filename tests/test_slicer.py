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