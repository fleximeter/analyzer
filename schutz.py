"""
File: schutz.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
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
import chart
import time
from decimal import Decimal
import os
import pathlib


def analyze():
    """
    Runs a salami slice analysis
    """
    # Path names
    path = "D:\\schutz\\"
    xml = os.path.join(path, "xml\\Schutz - Ein Kind Ist Uns Geboren - Full score - 01 Ein Kind Ist Uns Geboren.xml")
    voices = ["Soprano 1", "Soprano 2", "Alto", "Tenor 1", "Tenor 2", "Bass"]

    output = os.path.join(path, "register_analysis_files\\entire_piece.xlsx")
    output_general = os.path.join(path, "register_analysis_files\\statistics.xlsx")
    results_path = os.path.join(path, "register_analysis_files\\data.json")
    pathlib.Path(os.path.join(path, "register_analysis_files\\graphs")).mkdir(parents=True, exist_ok=True)

    # Record starting time
    start = time.time()
    use_cache = False

    # Analyze
    print("Analyzing entire piece...")
    results = None

    if use_cache:
        results = salami_slice_analyze.read_analysis_from_file(results_path)
    else:
        results = salami_slice_analyze.analyze(xml)
        salami_slice_analyze.write_analysis_to_file(results, results_path)

    salami_slice_analyze.write_general_report("Full piece", output_general, "w", results[0], results[0].lower_bound,
                                   results[0].upper_bound)
    salami_slice_analyze.write_report(output, results[0])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\chord_spacing_contours.xlsx"), ["chord_spacing_contour", "frequency", "duration"],
                               [results[0].chord_spacing_contour_frequency, results[0].chord_spacing_contour_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\psets.xlsx"), ["pset", "frequency", "duration"],
                               [results[0].pset_frequency, results[0].pset_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\pscs.xlsx"), ["psc", "frequency", "duration"],
                               [results[0].psc_frequency, results[0].psc_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\pcscs.xlsx"), ["sc", "frequency", "duration"],
                               [results[0].pcsc_frequency, results[0].pcsc_duration])

    # Make charts
    make_charts_general(results[0], path, voices)

    # Print elapsed time
    finish = time.time() - start
    print(f"\nTotal elapsed time: {int(finish / 60)} minutes, {round(finish % 60, 3)} seconds")


def make_charts_general(results, path, voices):
    """
    Makes general charts
    :param results: A Results object
    :param path: The file path
    :param voices: A list of voices
    :return:
    """
    chart.chart_cardinality(results, False, "Chord Cardinality Graph for Elliott Carter\u2019s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_m"))
    chart.chart_cardinality(results, True, "Chord Cardinality Graph for Elliott Carter\u2019s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_t"))
    chart.chart_cardinality_with_duplicates(results, False, "Chord Cardinality Graph for Elliott Carter\u2019s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_dup_m"))
    chart.chart_cardinality_with_duplicates(results, True, "Chord Cardinality Graph for Elliott Carter\u2019s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_dup_t"))
    chart.chart_pitch_onset(results, False, "Pitch Onsets in Elliott Carter\u2019s Fifth String Quartet", (6.5, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\onset_measure"))
    chart.chart_chord_spacing_index(results, False, "Chord Spacing Indices in Elliott Carter\u2019s Fifth String Quartet",
                                   (6.5, 3), os.path.join(path, "register_analysis_files\\graphs\\chord_spacing_index_m"))
    chart.chart_chord_spacing_index(results, True, "Chord Spacing Indices in Elliott Carter\u2019s Fifth String Quartet",
                                   (6.5, 3), os.path.join(path, "register_analysis_files\\graphs\\chord_spacing_index_t"))
    for i in range(len(voices)):
        chart.chart_pitch_onset(results, False, f"Pitch Onsets in Elliott Carter\u2019s Fifth String Quartet "
                                                   f"({voices[i]})", (6.5, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\onset_measure_{voices[i]}"), i)
    chart.chart_pitch_onset(results, True, "Pitch Onsets in Elliott Carter\u2019s Fifth String Quartet", (6.5, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\onset_time"))
    for i in range(len(voices)):
        chart.chart_pitch_onset(results, True, f"Pitch Onsets in Elliott Carter\u2019s Fifth String Quartet "
                                                  f"({voices[i]})", (6.5, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\onset_time_{voices[i]}"), i)
    chart.chart_pitch_duration(results, "Pitch Duration in Elliott Carter\u2019s Fifth String Quartet", (6.5, 3),
                               os.path.join(path, "register_analysis_files\\graphs\\pitch_duration"))
    for i in range(len(voices)):
        chart.chart_pitch_duration(results, f"Pitch Duration in Elliott Carter\u2019s Fifth String Quartet "
                                               f"({voices[i]})", (6.5, 3),
                                   os.path.join(path, f"register_analysis_files\\graphs\\pitch_duration_{voices[i]}"), i)
    chart.chart_pc_duration(results, "Pitch-Class Duration in Elliott Carter\u2019s Fifth String Quartet", (3, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\pc_duration"))
    for i in range(len(voices)):
        chart.chart_pc_duration(results, f"Pitch-Class Duration in Elliott Carter\u2019s Fifth String Quartet "
                                            f"({voices[i]})", (3, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\pc_duration_{voices[i]}"), i)


if __name__ == "__main__":
    print("################### Salami Slice Analyzer ####################\n" + \
          "Copyright (c) 2024 by Jeffrey Martin. All rights reserved.\nhttps://www.jeffreymartincomposer.com\n")
    analyze()
