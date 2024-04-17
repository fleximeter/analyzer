"""
File: carter.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functionality for analyzing register for Carter's fifth string quartet.
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
from fractions import Fraction
from decimal import Decimal
import os


def c_analyze():
    """
    Analyzes Carter's fifth string quartet without analyzing each section separately
    """
    xml = r"D:\carter_paper\xml\Carter String Quartet 5 - Full score - 01 Introduction.xml"
    output = r"D:\carter_paper\register_analysis_files\results_carter.xlsx"
    start = time.time()
    salami_slice_analyze.analyze(xml, output)
    finish = time.time() - start
    print(int(finish / 60), "minutes,", round(finish % 60, 3), "seconds")


def c_analyze_with_sections():
    """
    Analyzes Carter's fifth string quartet, as well as analyzing each section separately
    """
    # Sections
    measure_nos = [
        (1, 24), (25, 64), (65, 85), (86, 110), (111, 132), (133, 164), (165, 192), (193, 222), (223, 250),
        (251, 281), (282, 308), (309, 331)
    ]
    sections = [
        (1, 24), (25, 64), (65, 85), (86, 110), (111, 132), (133, 164), (165, 192), (193, 222), (223, 250),
        (251, 281), (282, 308), (309, 331), (1, 24), (25, 64), (65, 85), (86, 110), (111, 132), (133, 164),
        (165, 192), (193, 222), (223, 250), (251, 281), (282, 308), (309, 331)
    ]
    section_names = [
        "Introduction", "Giocoso", "Interlude I", "Lento espressivo", "Interlude II", "Presto scorrevole",
        "Interlude III", "Allegro energico", "Interlude IV", "Adagio sereno", "Interlude V", "Capriccioso"
    ]
    bound_prefs = [
        False, False, False, False, False, False, False, False, False, False, False, False,
        True, True, True, True, True, True, True, True, True, True, True, True
    ]
    voices = ["Violin 1", "Violin 2", "Viola", "Cello"]

    # Path names
    path = "D:\\carter_paper\\"
    xml = os.path.join(path, "xml\\Carter String Quartet 5 - Full score - 01 Introduction.xml")
    output = os.path.join(path, "register_analysis_files\\entire_piece.xlsx")
    output_general = os.path.join(path, "register_analysis_files\\statistics.xlsx")
    results_path = os.path.join(path, "register_analysis_files\\data.json")

    # the global output data
    output_global = []
    for i in range(12):
        cname = section_names[i].split(" ")
        c_path = os.path.join(path, f"register_analysis_files\\{i + 1}_")
        for j in range(len(cname) - 1):
            c_path += f"{cname[j]}_"
        c_path += f"{cname[len(cname) - 1]}_broad.xlsx"
        output_global.append(c_path)

    # the local output data
    output_local = []
    for i in range(12):
        cname = section_names[i].split(" ")
        c_path = os.path.join(path, f"register_analysis_files\\{i + 1}_")
        for j in range(len(cname) - 1):
            c_path += f"{cname[j]}_"
        c_path += f"{cname[len(cname) - 1]}_local.xlsx"
        output_local.append(c_path)

    # Record starting time
    start = time.time()
    use_cache = False

    # Analyze
    print("Analyzing entire piece...")
    results = None

    if use_cache:
        results = salami_slice_analyze.read_analysis_from_file(results_path)
    else:
        results = salami_slice_analyze.analyze_with_sections(xml, sections, bound_prefs)
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
    for i in range(1, len(output_global) + 1):
        salami_slice_analyze.write_general_report(f"Section {i} global", output_general, "a", results[i],
                                       results[0].lower_bound, results[0].upper_bound)
        salami_slice_analyze.write_report(output_global[i - 1], results[i])
        salami_slice_analyze.write_statistics(os.path.join(path,f"register_analysis_files\\chord_spacing_contours_{i}.xlsx"), ["chord_spacing_contour", "frequency", "duration"],
                                   [results[i].chord_spacing_contour_frequency, results[i].chord_spacing_contour_duration])
        salami_slice_analyze.write_statistics(os.path.join(path, f"register_analysis_files\\psets_{i}.xlsx"), ["pset", "frequency", "duration"],
                                   [results[i].pset_frequency, results[i].pset_duration])
        salami_slice_analyze.write_statistics(os.path.join(path, f"register_analysis_files\\pscs_{i}.xlsx"), ["psc", "frequency", "duration"],
                                   [results[i].psc_frequency, results[i].psc_duration])
        salami_slice_analyze.write_statistics(os.path.join(path, f"register_analysis_files\\pcscs_{i}.xlsx"), ["sc", "frequency", "duration"],
                                   [results[i].pcsc_frequency, results[i].pcsc_duration])

    for i in range(13, len(output_local) + 13):
        salami_slice_analyze.write_general_report(f"Section {i - 12} local", output_general, "a", results[i],
                                       results[0].lower_bound, results[0].upper_bound)
        salami_slice_analyze.write_report(output_local[i - 13], results[i])

    # Make charts
    make_charts_general(results[0], path, voices)

    for i in range(1, 13):
        make_charts_sections(results, i, path, voices, section_names)

    # for i in range(results[0].lower_bound, results[0].upper_bound + 1):
    #     total = Decimal(0)
    #     parts = Decimal(0)
    #     if i in results[0].pitch_duration:
    #         total = results[0].pitch_duration[i]
    #     for j in range(4):
    #         if i in results[0].pitch_duration_voices[j]:
    #             parts += results[0].pitch_duration_voices[j][i]
    #     print(f"{i}: {total - parts}")

    # Print elapsed time
    finish = time.time() - start
    print(f"\nTotal elapsed time: {int(finish / 60)} minutes, {round(finish % 60, 3)} seconds")


def c_analyze_reduction():
    """
    Analyzes a separate reduction for Section 12
    :return: None
    """
    path = "D:\\carter_paper\\"
    xml = os.path.join(path, "xml\\Section 12 Reduction - Full score - 01 12. Capriccioso.xml")
    output = os.path.join(path, "register_analysis_files\\reduction_section_12\\sec12_reduction.xlsx")
    output_general = os.path.join(path, "register_analysis_files\\reduction_section_12\\sec12_reduction_statistics.xlsx")
    results_path = os.path.join(path, "register_analysis_files\\data12.json")

    # Record starting time
    start = time.time()
    use_cache = False

    # Analyze
    print("Analyzing section 12 reduction...")
    results = None

    if use_cache:
        results = salami_slice_analyze.read_analysis_from_file(results_path)
    else:
        results = salami_slice_analyze.analyze(xml)
        salami_slice_analyze.write_analysis_to_file(results, results_path)

    salami_slice_analyze.write_general_report("Section 12", output_general, "w", results[0], results[0].lower_bound,
                                   results[0].upper_bound)
    salami_slice_analyze.write_report(output, results[0])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\reduction_section_12\\chord_spacing_contours.xlsx"),
                               ["chord_spacing_contour", "frequency", "duration"], [results[0].chord_spacing_contour_frequency, results[0].chord_spacing_contour_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\reduction_section_12\\psets.xlsx"),
                               ["pset", "frequency", "duration"], [results[0].pset_frequency, results[0].pset_duration])
    salami_slice_analyze.write_statistics(os.path.join(path, "register_analysis_files\\reduction_section_12\\pscs.xlsx"),
                               ["psc", "frequency", "duration"], [results[0].psc_frequency, results[0].psc_duration])

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
    chart.chart_cardinality(results, False, "Chord Cardinality Graph for Elliott Carter’s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_m"))
    chart.chart_cardinality(results, True, "Chord Cardinality Graph for Elliott Carter’s Fifth String Quartet",
                            size=(6.5, 3), path=os.path.join(path, "register_analysis_files\\graphs\\card_t"))
    chart.chart_pitch_onset(results, False, "Pitch Onsets in Elliott Carter’s Fifth String Quartet", (6.5, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\onset_measure"))
    chart.chart_chord_spacing_index(results, False, "Chord Spacing Indices in Elliott Carter's Fifth String Quartet",
                                   (6.5, 3), os.path.join(path, "register_analysis_files\\graphs\\chord_spacing_index_m"))
    chart.chart_chord_spacing_index(results, True, "Chord Spacing Indices in Elliott Carter's Fifth String Quartet",
                                   (6.5, 3), os.path.join(path, "register_analysis_files\\graphs\\chord_spacing_index_t"))
    for i in range(len(voices)):
        chart.chart_pitch_onset(results, False, f"Pitch Onsets in Elliott Carter’s Fifth String Quartet "
                                                   f"({voices[i]})", (6.5, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\onset_measure_{voices[i]}"), i)
    chart.chart_pitch_onset(results, True, "Pitch Onsets in Elliott Carter’s Fifth String Quartet", (6.5, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\onset_time"))
    for i in range(len(voices)):
        chart.chart_pitch_onset(results, True, f"Pitch Onsets in Elliott Carter’s Fifth String Quartet "
                                                  f"({voices[i]})", (6.5, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\onset_time_{voices[i]}"), i)
    chart.chart_pitch_duration(results, "Pitch Duration in Elliott Carter’s Fifth String Quartet", (6.5, 3),
                               os.path.join(path, "register_analysis_files\\graphs\\pitch_duration"))
    for i in range(len(voices)):
        chart.chart_pitch_duration(results, f"Pitch Duration in Elliott Carter’s Fifth String Quartet "
                                               f"({voices[i]})", (6.5, 3),
                                   os.path.join(path, f"register_analysis_files\\graphs\\pitch_duration_{voices[i]}"), i)
    chart.chart_pc_duration(results, "Pitch-Class Duration in Elliott Carter’s Fifth String Quartet", (3, 3),
                            os.path.join(path, "register_analysis_files\\graphs\\pc_duration"))
    for i in range(len(voices)):
        chart.chart_pc_duration(results, f"Pitch-Class Duration in Elliott Carter’s Fifth String Quartet "
                                            f"({voices[i]})", (3, 3),
                                os.path.join(path, f"register_analysis_files\\graphs\\pc_duration_{voices[i]}"), i)


def make_charts_sections(results, i, path, voices, section_names):
    """
    Makes charts
    :param results: A Results object
    :param i: The index of the Result
    :param path: The path of the register analysis files
    :param voices: A list of voices
    :param section_names: A list of section names
    :return:
    """
    # Create file names
    cname = section_names[i - 1].split(" ")
    cm_path = os.path.join(path, f"register_analysis_files\\graphs\\card_m_{i}_")
    ct_path = os.path.join(path, f"register_analysis_files\\graphs\\card_t_{i}_")
    csim_path = os.path.join(path, f"register_analysis_files\\graphs\\csi_m_{i}_")
    csit_path = os.path.join(path, f"register_analysis_files\\graphs\\csi_t_{i}_")
    om_path = os.path.join(path, f"register_analysis_files\\graphs\\onset_m_{i}_")
    ot_path = os.path.join(path, f"register_analysis_files\\graphs\\onset_t_{i}_")
    dp_path = os.path.join(path, f"register_analysis_files\\graphs\\dur_pitch_{i}_")
    dpc_path = os.path.join(path, f"register_analysis_files\\graphs\\dur_pc_{i}_")
    for j in range(len(cname) - 1):
        cm_path += f"{cname[j]}_"
        ct_path += f"{cname[j]}_"
        om_path += f"{cname[j]}_"
        ot_path += f"{cname[j]}_"
        dp_path += f"{cname[j]}_"
        dpc_path += f"{cname[j]}_"
    cm_path += cname[len(cname) - 1]
    ct_path += cname[len(cname) - 1]
    om_path += cname[len(cname) - 1]
    ot_path += cname[len(cname) - 1]
    dp_path += cname[len(cname) - 1]
    dpc_path += cname[len(cname) - 1]

    # Create charts
    chart.chart_cardinality(results[i], False, f"Chord Cardinality Graph for Section {i} – {section_names[i - 1]}",
                            path=cm_path)
    chart.chart_cardinality(results[i], True, f"Chord Cardinality Graph for Section {i} – {section_names[i - 1]}",
                            path=ct_path)
    chart.chart_chord_spacing_index(results[i], False, f"Chord Spacing Indices in Section {i} – {section_names[i - 1]}",
                                   path=csim_path)
    chart.chart_chord_spacing_index(results[i], True, f"Chord Spacing Indices in Section {i} – {section_names[i - 1]}",
                                   path=csit_path)
    chart.chart_pitch_onset(results[i], False, f"Pitch Onsets in Section {i} – {section_names[i - 1]}",
                            path=om_path)
    for j in range(len(voices)):
        chart.chart_pitch_onset(results[i], False, f"Pitch Onsets in Section {i} ({voices[j]}) – "
                                                   f"{section_names[i - 1]}", path=f"{om_path}_{voices[j]}", voice=j)
    chart.chart_pitch_onset(results[i], True, f"Pitch Onsets in Section {i} – {section_names[i - 1]}",
                            path=ot_path)
    for j in range(len(voices)):
        chart.chart_pitch_onset(results[i], True, f"Pitch Onsets in Section {i} ({voices[j]}) – {section_names[i - 1]}",
                                path=f"{ot_path}_{voices[j]}", voice=j)
    chart.chart_pitch_duration(results[i], f"Pitch Durations in Section {i} – {section_names[i - 1]}",
                               path=dp_path)
    for j in range(len(voices)):
        chart.chart_pitch_duration(results[i], f"Pitch Durations in Section {i} ({voices[j]}) – {section_names[i - 1]}",
                                   path=f"{dp_path}_{voices[j]}", voice=j)
    chart.chart_pc_duration(results[i], f"Pitch-Class Durations in Section {i} – {section_names[i - 1]}",
                            path=dpc_path)
    for j in range(len(voices)):
        chart.chart_pc_duration(results[i], f"Pitch-Class Durations in Section {i} ({voices[j]}) – "
                                            f"{section_names[i - 1]}", path=f"{dpc_path}_{voices[j]}", voice=j)


def metric_modulation():
    # Carter Quartet 5
    m = {}

    # q = 72
    m[24] = Fraction(Fraction(3, 4), 1)  # q = 96
    m[25] = Fraction(Fraction(1, 4), Fraction(1, 6))  # q = 64
    m[45] = Fraction(Fraction(1, 8), Fraction(1, 7))  # q = 73+ (512/7)
    m[65] = Fraction(Fraction(4, 3), 1)  # q = 55- (384/7)
    m[71] = Fraction(Fraction(4, 7), 1)  # q = 96
    m[77] = Fraction(Fraction(8, 5), 1)  # q = 60
    m[123] = Fraction(Fraction(5, 4), 1)  # q = 48
    m[127] = Fraction(Fraction(Fraction(5, 2), 3), 2)  # h = 57.6
    m[174] = Fraction(Fraction(8, 5), 1)  # q = 72
    m[231] = Fraction(Fraction(2, 3), 1)  # q = 108
    m[238] = Fraction(Fraction(3, 2), 1)  # q = 72
    m[241] = Fraction(Fraction(3, 4), Fraction(1, 2))  # q = 48
    m[282] = Fraction(Fraction(1, 2), 1)  # q = 96
    m[308] = Fraction(Fraction(4, 5), Fraction(1, 2))  # q = 60

    for i in m:
        print(i, m[i])


if __name__ == "__main__":
    print("################### Salami Slice Analyzer ####################\n" + \
          "Copyright (c) 2024 by Jeffrey Martin. All rights reserved.\nhttps://www.jeffreymartincomposer.com\n")
    c_analyze_with_sections()
    # c_analyze_reduction()
    # metric_modulation()
