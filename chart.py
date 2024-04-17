"""
File: chart.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functionality for making analysis charts.
Copyright (c) 2022 by Jeff Martin.

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

import matplotlib.pyplot
import numpy


def chart_cardinality(results, x_axis_time=False, title="Cardinality Chart", size=(8, 6), path=None):
    """
    Makes a step plot of cardinality
    :param results: A Results object
    :param x_axis_time: Whether or not to use time as the x axis (default is measure number)
    :param title: The title of the plot
    :param size: The size of the plot
    :param path: A path to save the chart
    :return: None
    """
    matplotlib.pyplot.clf()
    matplotlib.pyplot.rcParams["font.family"] = "Academico"
    ps_s = []
    x = []
    position_time = results.start_time
    for s in results.slices:
        if x_axis_time:
            x.append(position_time / 60)
            position_time += s.duration
        else:
            position = s.measure + float(s.start_position / s.time_signature.barDuration.quarterLength)
            x.append(position)
        if s.p_cardinality > 0:
            ps_s.append(s.p_cardinality)
        else:
            ps_s.append(0)
    fig = matplotlib.pyplot.figure(figsize=size, dpi=600)
    matplotlib.pyplot.rc("font", size=8)
    matplotlib.pyplot.rc("axes", titlesize=8)
    matplotlib.pyplot.rc("axes", labelsize=8)
    matplotlib.pyplot.rc("axes", linewidth=0.5)
    matplotlib.pyplot.rc("xtick", labelsize=8)
    matplotlib.pyplot.rc("xtick.major", width=0.5)
    matplotlib.pyplot.rc("ytick", labelsize=8)
    matplotlib.pyplot.rc("ytick.major", width=0.5)
    ax = fig.add_subplot(111)
    ax.step(x, ps_s, color="#555555", linewidth=0.3, markersize=1)
    # matplotlib.pyplot.title(title, fontsize=18)
    if x_axis_time:
        matplotlib.pyplot.xlabel("Time (minutes)")
    else:
        matplotlib.pyplot.xlabel("Measure No.")
    matplotlib.pyplot.ylabel("Chord Cardinality")
    matplotlib.pyplot.tight_layout()
    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()


def chart_chord_spacing_index(results, x_axis_time=False, title="Chord Spacing Index Chart", size=(8, 6), path=None):
    """
    Makes a step plot of chord spacing index
    :param results: A Results object
    :param x_axis_time: Whether or not to use time as the x axis (default is measure number)
    :param title: The title of the plot
    :param size: The size of the plot
    :param path: A path to save the chart
    :return: None
    """
    matplotlib.pyplot.clf()
    matplotlib.pyplot.rcParams["font.family"] = "Academico"
    ps_s = []
    x = []
    position_time = results.start_time
    for s in results.slices:
        if x_axis_time:
            x.append(position_time / 60)
            position_time += s.duration
        else:
            position = s.measure + float(s.start_position / s.time_signature.barDuration.quarterLength)
            x.append(position)
        if s.p_cardinality > 0:
            ps_s.append(s.chord_spacing_index)
        else:
            ps_s.append(numpy.nan)
    fig = matplotlib.pyplot.figure(figsize=size, dpi=600)
    matplotlib.pyplot.rc("font", size=8)
    matplotlib.pyplot.rc("axes", titlesize=8)
    matplotlib.pyplot.rc("axes", labelsize=8)
    matplotlib.pyplot.rc("axes", linewidth=0.5)
    matplotlib.pyplot.rc("xtick", labelsize=8)
    matplotlib.pyplot.rc("xtick.major", width=0.5)
    matplotlib.pyplot.rc("ytick", labelsize=8)
    matplotlib.pyplot.rc("ytick.major", width=0.5)
    ax = fig.add_subplot(111)
    ax.step(x, ps_s, color="#555555", linewidth=0.3, markersize=1)
    # matplotlib.pyplot.title(title, fontsize=18)
    if x_axis_time:
        matplotlib.pyplot.xlabel("Time (minutes)")
    else:
        matplotlib.pyplot.xlabel("Measure No.")
    matplotlib.pyplot.ylabel("Chord Spacing Index")
    matplotlib.pyplot.tight_layout()
    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()


def chart_pitch_onset(results, x_axis_time=False, title="Pitch Onset Graph", size=(8, 6), path=None, voice=None):
    """
    Makes a pitch onset chart
    :param results: A Results object
    :param x_axis_time: Whether or not to use time as the x axis (default is measure number)
    :param title: The title of the plot
    :param size: The size of the plot
    :param path: A path to save the chart
    :param voice: The number of the voice (if None, charts all voices)
    :return: None
    """
    matplotlib.pyplot.clf()
    matplotlib.pyplot.rcParams["font.family"] = "Academico"
    pitches = [[] for i in range(results.max_p_count)]
    x = []
    position_time = results.start_time
    for s in range(len(results.slices)):
        if x_axis_time:
            x.append(position_time / 60)
            position_time += results.slices[s].duration
        else:
            position = results.slices[s].measure
            position += float(results.slices[s].start_position / results.slices[s].time_signature.barDuration.quarterLength)
            x.append(position)

        # Add the pitches in the pseg of the current slice to the plot
        for i in range(len(pitches)):
            if voice is None:
                if i < len(results.slices[s].pseg):
                    if s - 1 >= 0:
                        if results.slices[s].pseg[i] not in results.slices[s - 1].pset:
                            pitches[i].append(results.slices[s].pseg[i].p)
                        else:
                            pitches[i].append(numpy.nan)
                    else:
                        pitches[i].append(results.slices[s].pseg[i].p)
                else:
                    pitches[i].append(numpy.nan)
            else:
                if voice < len(results.slices[s].psegs):
                    if i < len(results.slices[s].psegs[voice]):
                        if s - 1 >= 0 and voice < len(results.slices[s - 1].psets):
                            if results.slices[s].psegs[voice][i] not in results.slices[s - 1].psets[voice]:
                                pitches[i].append(results.slices[s].psegs[voice][i].p)
                            else:
                                pitches[i].append(numpy.nan)
                        else:
                            pitches[i].append(results.slices[s].psegs[voice][i].p)
                    else:
                        pitches[i].append(numpy.nan)
                else:
                    pitches[i].append(numpy.nan)
    draw = matplotlib.pyplot.figure(figsize=size, dpi=600)
    matplotlib.pyplot.rc("font", size=8)
    matplotlib.pyplot.rc("axes", titlesize=8)
    matplotlib.pyplot.rc("axes", labelsize=8)
    matplotlib.pyplot.rc("axes", linewidth=0.5)
    matplotlib.pyplot.rc("xtick", labelsize=8)
    matplotlib.pyplot.rc("xtick.major", width=0.5)
    matplotlib.pyplot.rc("ytick", labelsize=8)
    matplotlib.pyplot.rc("ytick.major", width=0.5)
    axes = draw.add_subplot(111)
    for i in range(len(pitches)):
        axes.scatter(x, pitches[i], s=2, c='#222222', marker='.')
    ytick_tick = []
    ytick_labels = []
    for i in range(results.lower_bound, results.upper_bound + 1):
        if i % 12 == 0:
            ytick_tick.append(i)
            ytick_labels.append(f"C{i // 12 + 4}/{i}")
    matplotlib.pyplot.yticks(ytick_tick, ytick_labels)
    # matplotlib.pyplot.title(title, fontsize=18)
    if x_axis_time:
        matplotlib.pyplot.xlabel("Time (minutes)")
    else:
        matplotlib.pyplot.xlabel("Measure No.")
    matplotlib.pyplot.ylabel("Pitch")
    matplotlib.pyplot.tight_layout()
    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()


def chart_pitch_duration(results, title="Pitch Duration Graph", size=(8, 6), path=None, voice=None):
    """
    Makes a pitch duration bar chart
    :param results: A Results object
    :param title: The title of the chart
    :param size: The size of the plot
    :param path: A path to save the chart
    :param voice: The number of the voice (if None, charts all voices)
    :return: None
    """
    matplotlib.pyplot.clf()
    matplotlib.pyplot.rcParams["font.family"] = "Academico"
    x = [i for i in range(results.pitch_lowest, results.pitch_highest + 1)]
    y = []
    for p in x:
        if voice is None:
            if p in results.pitch_duration:
                y.append(float(results.pitch_duration[p]))
            else:
                y.append(0)
        else:
            if p in results.pitch_duration_voices[voice]:
                y.append(float(results.pitch_duration_voices[voice][p]))
            else:
                y.append(0)
    draw = matplotlib.pyplot.figure(figsize=size, dpi=600)
    axes = draw.add_subplot(111)
    axes.bar(x, y)
    xtick_tick = []
    xtick_labels = []
    for i in range(results.lower_bound, results.upper_bound + 1):
        if i % 12 == 0:
            xtick_tick.append(i)
            xtick_labels.append(f"C{i // 12 + 4}/{i}")
    matplotlib.pyplot.xticks(xtick_tick, xtick_labels)
    matplotlib.pyplot.title(title, fontsize=18)
    matplotlib.pyplot.xlabel("Pitch")
    matplotlib.pyplot.ylabel("Duration (seconds)")

    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()


def chart_pc_duration(results, title="Pitch-Class Duration Graph", size=(8, 6), path=None, voice=None):
    """
    Makes a pitch-class duration bar chart
    :param results: A Results object
    :param title: The title of the chart
    :param size: The size of the plot
    :param path: A path to save the chart
    :param voice: The number of the voice (if None, charts all voices)
    :return: None
    """
    matplotlib.pyplot.clf()
    matplotlib.pyplot.rcParams["font.family"] = "Academico"
    x = [i for i in range(12)]
    y = []
    for pc in x:
        if voice is None:
            if pc in results.pc_duration:
                y.append(float(results.pc_duration[pc]))
            else:
                y.append(0)
        else:
            if pc in results.pc_duration_voices[voice]:
                y.append(float(results.pc_duration_voices[voice][pc]))
            else:
                y.append(0)
    draw = matplotlib.pyplot.figure(figsize=size, dpi=600)
    axes = draw.add_subplot(111)
    axes.bar(x, y)
    xtick_tick = [i for i in range(12)]
    xtick_labels = [
        "C\n0",
        "C♯/D♭\n1",
        "D\n2",
        "D♯/E♭\n3",
        "E\n4",
        "F\n5",
        "F♯/G♭\n6",
        "G\n7",
        "G♯/A♭\n8",
        "A\n9",
        "A♯/B♭\n10",
        "B\n11"
    ]
    matplotlib.pyplot.xticks(xtick_tick, xtick_labels)
    matplotlib.pyplot.title(title, fontsize=18)
    matplotlib.pyplot.xlabel("Pitch-Class")
    matplotlib.pyplot.ylabel("Duration (seconds)")

    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()


def chart_spaces(results, size=(8, 6), path=None):
    """
    Makes a pitch space stacked area chart
    :param results: A Results object
    :param size: The size of the plot
    :param path: A path to save the chart
    :return: None
    """
    matplotlib.pyplot.clf()
    x = []
    ins = []
    lns = []
    ns = []
    ps = []
    uns = []
    for s in results.slices:
        position = s.measure + float(s.start_position / s.time_signature.barDuration.quarterLength)
        x.append(position)
        if s.ins is not None:
            ins.append(s.ins)
        else:
            ins.append(0)
        if s.ps is not None:
            ps.append(s.ps)
        else:
            ps.append(0)
        if s.uns is None:
            lns.append(0)
            uns.append(0)
            ns.append(results.lps_card)
        else:
            lns.append(s.lns)
            uns.append(s.uns)
            ns.append(0)
    """
    matplotlib.pyplot.step(x, uns, fillstyle="full")
    matplotlib.pyplot.step(x, ps)
    """
    labels = ["NS", "LNS", "INS", "PS", "UNS"]
    colors = ["#000000", "#333333", "#555555", "#DDDDDD", "#333333"]
    fig, ax = matplotlib.pyplot.subplots()
    matplotlib.pyplot.figure(figsize=size, dpi=600)
    ax.stackplot(x, ns, lns, ins, ps, uns, labels=labels, colors=colors)
    ax.legend(loc="upper left")
    if path is None:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(path)
        matplotlib.pyplot.close()
