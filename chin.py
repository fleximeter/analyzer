"""
File: chin.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functionality for analyzing register for Chin's Piano Etude (In C).
Copyright (c) 2024 by Jeff Martin.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import salami_slice_analyze
import rhythm_analyze
import chart
import time
import os
from fractions import Fraction


def c_analyze():
    """
    Analyzes Chin's "In C" without analyzing each section separately
    """
    # A map of tempos for the piece
    tempo_map = {1: 84, 26: 108}

    # Path names
    path = "D:\\chin_paper\\"
    path_laptop = "C:\\Users\\jeffr\\chin_paper\\"
    # path = path_laptop
    xml = os.path.join(path, "chin_etude_1_6staff.musicxml")
    xml_no_grace_notes = os.path.join(path, "chin_etude_1_6staff_no_grace_notes.musicxml")
    output = os.path.join(path, "analysis\\entire_piece.xlsx")
    output_general = os.path.join(path, "analysis\\statistics.xlsx")
    results_path = os.path.join(path, "analysis\\data.json")
    
    # Record starting time
    start = time.time()

    # Analyze
    print("Analyzing entire piece...")
    results = None
    
    results = salami_slice_analyze.analyze(xml_no_grace_notes, tempo_map=tempo_map)
    results_staff = [salami_slice_analyze.analyze(xml_no_grace_notes, staff_indices=[i], tempo_map=tempo_map) for i in range(0, 6)]

    salami_slice_analyze.write_general_report("Full piece", output_general, "w", results[0], results[0].lower_bound,
                                   results[0].upper_bound)
    salami_slice_analyze.write_report(output, results[0])
    salami_slice_analyze.write_statistics(os.path.join(path, "analysis", "chord_spacing_contours.xlsx"), ["chord_spacing_contour", "frequency", "duration"],
                               [results[0].chord_spacing_contour_frequency, results[0].chord_spacing_contour_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "analysis", "psets.xlsx"), ["pset", "frequency", "duration"],
                               [results[0].pset_frequency, results[0].pset_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "analysis", "pscs.xlsx"), ["psc", "frequency", "duration"],
                               [results[0].psc_frequency, results[0].psc_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "analysis", "pcscs.xlsx"), ["sc", "frequency", "duration"],
                               [results[0].pcsc_frequency, results[0].pcsc_duration])

    make_charts_general(results[0], os.path.join(path, "analysis", "graphs"))
    
    for i, result in enumerate(results_staff):
        salami_slice_analyze.write_report(os.path.join(path, "analysis", f"entire_piece_staff{i+1}.xlsx"), result[0])
        salami_slice_analyze.write_statistics(os.path.join(path, "analysis", f"psets_staff{i+1}.xlsx"), ["pset", "frequency", "duration"],
                                [result[0].pset_frequency, result[0].pset_duration])
        make_charts_specific(result[0], os.path.join(path, "analysis", f"graphs_staff{i+1}"))
    
    # Perform IOI analysis
    ioi_analysis = rhythm_analyze.part_ioi_analyzer(xml_no_grace_notes)
    chart.chart_voice_ioi(ioi_analysis, [f"Voice {i+1}" for i in range(len(ioi_analysis))], ("#0066ff", "#ff9900", "#00cc00", "#ff0000", "#cc66ff", "#cc9900"),
                          "IOI Voice Graph for Unsuk Chin\u2019s \u201cIn C\u201d", (10, 6), os.path.join(path, "analysis", "graphs", "ioi_voice_graph"))
    chart.chart_voice_ioi([ioi_analysis[1], ioi_analysis[3]], (f"Voice 2", "Voice 4"), ("#ff9900", "#ff0000"),
                          "IOI Voice Graph for Unsuk Chin\u2019s \u201cIn C\u201d (Voices 2, 4)", (10, 6), os.path.join(path, "analysis", "graphs", "ioi_voice_graph_24"))
    rhythm_analyze.write_ioi_analysis_to_file(os.path.join(path, "analysis", "ioi.xlsx"), ioi_analysis)
    
    # Find relationships
    print("2-1 Relationships")
    print(rhythm_analyze.find_rhythm_ratio_successsion(ioi_analysis, (2, 1)))
    print("2-1-2-1 Relationships")
    print(rhythm_analyze.find_rhythm_ratio_successsion(ioi_analysis, (2, 1, 2, 1)))
    print("2-1-2-1-2-1 Relationships")
    print(rhythm_analyze.find_rhythm_ratio_successsion(ioi_analysis, (2, 1, 2, 1, 2, 1)))
    print("2-1-2-1-2-1-2-1 Relationships")
    print(rhythm_analyze.find_rhythm_ratio_successsion(ioi_analysis, (2, 1, 2, 1, 2, 1, 2, 1)))
    print("e-s Relationships")
    print(rhythm_analyze.find_duration_successsion(ioi_analysis, (Fraction(1, 2), Fraction(1, 4))))
    print("e-s-e-s Relationships")
    print(rhythm_analyze.find_duration_successsion(ioi_analysis, (Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4))))
    print("e-s-e-s-e-s Relationships")
    print(rhythm_analyze.find_duration_successsion(ioi_analysis, (Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4))))
    print("e-s-e-s-e-s-e-s Relationships")
    print(rhythm_analyze.find_duration_successsion(ioi_analysis, (Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4), Fraction(1, 2), Fraction(1, 4))))
    
    # Print elapsed time
    finish = time.time() - start
    print(f"\nTotal elapsed time: {int(finish / 60)} minutes, {round(finish % 60, 3)} seconds")


def make_charts_general(results, path):
    """
    Makes general charts
    :param results: A Results object
    :param path: The file path
    :return:
    """
    chart.chart_cardinality(results, False, "Chord Cardinality Graph for Unsuk Chin\u2019s \u201cIn C\u201d",
                            size=(6.5, 3), path=os.path.join(path, "card_m"))
    chart.chart_cardinality(results, True, "Chord Cardinality Graph for Unsuk Chin\u2019s \u201cIn C\u201d",
                            size=(6.5, 3), path=os.path.join(path, "card_t"))
    chart.chart_pitch_onset(results, False, "Pitch Onsets in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                            os.path.join(path, "onset_measure"))
    chart.chart_chord_spacing_index(results, False, "Chord Spacing Indices in Unsuk Chin\u2019s \u201cIn C\u201d",
                                   (6.5, 3), os.path.join(path, "chord_spacing_index_m"))
    chart.chart_chord_spacing_index(results, True, "Chord Spacing Indices in Unsuk Chin\u2019s \u201cIn C\u201d",
                                   (6.5, 3), os.path.join(path, "chord_spacing_index_t"))
    chart.chart_pitch_onset(results, True, "Pitch Onsets in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                            os.path.join(path, "onset_time"))
    chart.chart_pitch_duration(results, "Pitch Duration in Unsuk Chin\u2019s \u201cIn C\u201d", (6, 3.5),
                               os.path.join(path, "pitch_duration"))
    chart.chart_pc_duration(results, "Pitch-Class Duration in Unsuk Chin\u2019s \u201cIn C\u201d", (6, 3.5),
                            os.path.join(path, "pc_duration"))
                            

def make_charts_specific(results, path):
    """
    Make specific charts
    :param results: A Results object
    :param path: The file path
    :return:
    """
    chart.chart_cardinality(results, False, "Pset Cardinality Graph for Unsuk Chin\u2019s \u201cIn C\u201d",
                            size=(6.5, 3), path=os.path.join(path, "card_m"))
    chart.chart_cardinality(results, True, "Pset Cardinality Graph for Unsuk Chin\u2019s \u201cIn C\u201d",
                            size=(6.5, 3), path=os.path.join(path, "card_t"))
    chart.chart_pitch_onset(results, False, "Pitch Onsets in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                            os.path.join(path, "onset_measure"))
    chart.chart_pitch_onset(results, True, "Pitch Onsets in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                            os.path.join(path, "onset_time"))
    chart.chart_pitch_duration(results, "Pitch Duration in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                               os.path.join(path, "pitch_duration"))
    chart.chart_pc_duration(results, "Pitch-Class Duration in Unsuk Chin\u2019s \u201cIn C\u201d", (6.5, 3),
                            os.path.join(path, "pc_duration"))
    

if __name__ == "__main__":
    print("################### Salami Slice Analyzer ####################\n" + \
          "Copyright (c) 2024 by Jeffrey Martin. All rights reserved.\nhttps://www.jeffreymartincomposer.com\n")
    c_analyze()
