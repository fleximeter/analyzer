"""
File: salami_slice_analyze.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functions for analyzing salami slices (distinct vertical 
simultaneities) of a score
Copyright (c) 2022, 2024 by Jeff Martin.

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


import json
import music21
import pandas as pd
from salami_slice import SalamiSlice
from results import Results
from fractions import Fraction
from pctheory import pitch, pcset
from decimal import Decimal


def analyze(input_xml, starting_measure_num=None, ending_measure_num=None, use_local=False, staff_indices=None, tempo_map=None):
    """
    Performs a vertical analysis on the given stream and writes a report to CSV
    :param input_xml: The musicxml file to analyze
    :param starting_measure_num: The first measure to analyze
    :param ending_measure_num: The last measure to analyze
    :param use_local: Whether or not to use local bounds for register analysis
    :param staff_indices: Whether to only analyze some of the staves or analyze 
    the whole score (a value of None means to analyze the whole score)
    :param tempo_map: A map of tempos to use to override the tempos in the score. If None, tempos will 
    be read from the score. The map should be structured with measures as keys, and tempos as values.
    :return: A Results object containing the results of the analysis
    """
    stream = music21.converter.parse(input_xml)
    parts = []
    
    if staff_indices is not None:
        i = 0
        for item in stream:
            if type(item) == music21.stream.Part or type(item) == music21.stream.PartStaff:
                if i in staff_indices:
                    parts.append(item)
                i += 1
    
    else:            
        for item in stream:
            if type(item) == music21.stream.Part or type(item) == music21.stream.PartStaff:
                parts.append(item)
    
    results = slice_parts(parts, get_slice_num(parts), [], [use_local], starting_measure_num, ending_measure_num, tempo_map)
    return results


def analyze_corpus(name, starting_measure_num=None, ending_measure_num=None, use_local=False):
    """
    Performs a vertical analysis on the given stream and writes a report to CSV
    :param name: The musicxml file in the music21 corpus to analyze
    :param starting_measure_num: The first measure to analyze
    :param ending_measure_num: The last measure to analyze
    :param use_local: Whether or not to use local bounds for register analysis
    :return: A Results object containing the results of the analysis
    """
    stream = music21.corpus.parse(name)
    parts = []
    for item in stream:
        if type(item) == music21.stream.Part or type(item) == music21.stream.PartStaff:
            parts.append(item)
    results = slice_parts(parts, get_slice_num(parts), [], [use_local], starting_measure_num, ending_measure_num)
    return results[0]


def analyze_with_sections(input_xml, section_divisions, use_local, starting_measure_num=None, ending_measure_num=None, tempo_map=None):
    """
    Performs a vertical analysis on the given stream and writes a report to CSV
    :param input_xml: The musicxml file to analyze
    :param section_divisions: A list of section divisions
    :param use_local: Whether or not to use local bounds for register analysis
    :param tempo_map: A map of tempos to use to override the tempos in the score. If None, tempos will 
    be read from the score. The map should be structured with measures as keys, and tempos as values.
    :return: A list of Results objects containing the results of the analysis.
    Index 0 is a complete analysis, and the remaining indices are section analyses
    in the order in which they were provided.
    """
    stream = music21.converter.parse(input_xml)
    parts = []
    for item in stream:
        if type(item) == music21.stream.Part or type(item) == music21.stream.PartStaff:
            parts.append(item)
    return slice_parts(parts, get_slice_num(parts), section_divisions, use_local, starting_measure_num, ending_measure_num, tempo_map=tempo_map)


def clean_slices(salami_slices, match_tempo=False, sections=None):
    """
    Cleans up a list of v_slices
    :param slices: A list of v_slices
    :param match_tempo: Whether or not to force tempo match
    :param sections: A list of section divisions
    """
    # Remove duplicate slices, and update durations
    i = 1
    while i < len(salami_slices):
        equal = True
        if match_tempo and salami_slices[i]._tempo != salami_slices[i - 1]._tempo:
            equal = False
        elif salami_slices[i].pitchseg != salami_slices[i - 1].pitchseg:
            equal = False
        elif sections is not None:
            if salami_slices[i].measure in sections and salami_slices[i - 1].measure < salami_slices[i].measure:
                equal = False
        if equal:
            salami_slices[i - 1].duration += salami_slices[i].duration
            salami_slices[i - 1].quarter_duration += salami_slices[i].quarter_duration
            del salami_slices[i]
        else:
            i += 1
    
    # Calculate IOI
    for j, salami_slice in enumerate(salami_slices):
        if salami_slice.pset_cardinality > 0:
            if j + 1 < len(salami_slices) and salami_slices[j+1].pset_cardinality == 0:
                salami_slice._ioi_in_seconds = salami_slices[j].duration + salami_slices[j+1].duration
            else:
                salami_slice._ioi_in_seconds = salami_slice.duration


def factor(n):
    """
    Factors a positive integer
    :param n: An integer
    :returns: A list of factors, in sorted order, including duplicates
    """
    factors = []
    d = 1
    while d <= int(n ** 0.5):
        if n % d == 0:
            factors.append(d)
            n //= d
        else:
            d += 1
        if d == 1:
            d += 1
        # if d > int(n ** 0.5):
        #     factors.append(n)
    factors.append(n)
    # factors.sort()
    return factors


#done
def get_bounds(slices):
    """
    Gets the upper and lower bounds of a list of v_slices
    :param slices: A list of v_slices
    :return: The lower and upper bounds as a tuple. The lower bound is index 0,
    and the upper bound is index 1. If the slices contain no pitches, each of
    the bounds will be None.
    """

    lower_bound = None
    upper_bound = None

    for i in range(0, len(slices)):
        if (lower_bound is None or upper_bound is None) and len(slices[i].pseg) > 0:
            lower_bound = slices[i].pseg[0].p
            upper_bound = slices[i].pseg[len(slices[i].pseg) - 1].p
        if len(slices[i].pseg) > 0 and lower_bound is not None and upper_bound is not None:
            if slices[i].pseg[0].p < lower_bound:
                lower_bound = slices[i].pseg[0].p
            if slices[i].pseg[len(slices[i].pseg) - 1].p > upper_bound:
                upper_bound = slices[i].pseg[len(slices[i].pseg) - 1].p

    return lower_bound, upper_bound


#done
def get_piece_bounds(parts):
    """
    Determines the lower and upper bounds of a piece
    :param parts: A list of parts
    :return: The lower and upper bounds as a tuple.
    """
    lower = None
    upper = None
    for part in parts:
        for item in part:
            if type(item) == music21.stream.Measure:
                for item2 in item:
                    if type(item2) == music21.stream.Voice:
                        for item3 in item2:
                            if type(item3) == music21.note.Note or type(item3) == music21.chord.Chord:
                                for pitch in item3.pitches:
                                    if lower is None:
                                        lower = pitch.midi - 60
                                    if upper is None:
                                        upper = pitch.midi - 60
                                    if lower > pitch.midi - 60:
                                        lower = pitch.midi - 60
                                    if upper < pitch.midi - 60:
                                        upper = pitch.midi - 60
                    elif type(item2) == music21.note.Note or type(item2) == music21.chord.Chord:
                        for pitch in item2.pitches:
                            if lower is None:
                                lower = pitch.midi - 60
                            if upper is None:
                                upper = pitch.midi - 60
                            if lower > pitch.midi - 60:
                                lower = pitch.midi - 60
                            if upper < pitch.midi - 60:
                                upper = pitch.midi - 60

    return lower, upper


#done
def get_slice_num(parts):
    """
    Determines the number of slices per quarter note based on subdivisions of the note.
    :param parts: A stream of parts
    :returns: The number of slices per quarter note
    """
    # A collection of all the unique denominators we find
    denominators = {}
    denominators_list = []

    # Find all the unique denominators
    for part in parts:
        for stream in part:
            if type(stream) == music21.stream.Measure:
                for item in stream:
                    if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(item) == \
                            music21.chord.Chord:
                        ql = item.duration.quarterLength
                        if type(item.duration.quarterLength) != Fraction:
                            ql = Fraction(item.duration.quarterLength)
                        if ql.denominator not in denominators:
                            denominators[ql.denominator] = True

    # Get the LCM and return it. This is the number of slices per quarter note that we need.
    for item in denominators:
        denominators_list.append(item)
    # print(lcm(denominators_list))
    return lcm(denominators_list)


#done
def lcm(integers):
    """
    Computes the LCM of a list of positive integers
    :param integers: A list of positive integers
    :return: The LCM
    """
    factors = {}  # A dictionary of individual factors and their multiplicities
    multiple = 1  # The LCM

    for num in integers:
        cur_factors = factor(num)  # The factors of the current number
        current = 1  # The current factor we are considering
        count = 0  # The number of occurrences of that factor
        for i in range(len(cur_factors)):
            # If we found another occurrence of that factor, increase the count
            if cur_factors[i] == current:
                count += 1
            # Otherwise record the count and move on
            else:
                if current not in factors:
                    factors[current] = count
                elif factors[current] < count:
                    factors[current] = count
                current = cur_factors[i]
                count = 1
            # If we are done, record the count of the last factor
            if i + 1 == len(cur_factors):
                if current not in factors:
                    factors[current] = count
                elif factors[current] < count:
                    factors[current] = count

    # Compute the LCM
    for item in factors:
        multiple *= item ** factors[item]
    # print(multiple)
    return multiple


#done
def set_slice_bounds(slices, bounds):
    """
    Sets the bounds of a list of v_slices
    :param slices: A list of v_slices
    :param bounds: A tuple with the lower and upper bounds
    """

    for i in range(len(slices)):
        slices[i].lower_bound = bounds[0]
        slices[i].upper_bound = bounds[1]


#done
def slice_parts(parts, n, section_divisions, use_local, starting_measure_num=None, ending_measure_num=None, tempo_map=None):
    """
    Takes n vertical slices of each beat from each of the parts. Note that beats are always quarter notes
    in music21. The parts do not need to have the same time signature for each measure: each slice is taken
    independently of the time signature. The parts do not even need to have the same number of total beats.
    However, it is assumed that a quarter note in any given part is equal in duration to a quarter note in
    any other part (this means that all parts must share the same tempo for a quarter note).
    :param parts: A list of parts
    :param n: The number of slices per quarter note
    :param section_divisions: A list of section divisions
    :param use_local: Whether or not to use local bounds for register analysis
    :param starting_measure_num: The first measure to analyze (None means start at the beginning)
    :param ending_measure_num: The last measure to analyze (None means analyze to the end)
    :param tempo_map: A map of tempos to use to override the tempos in the score. If None, tempos will 
    be read from the score. The map should be structured with measures as keys, and tempos as values.
    :return: A list of v_slices
    """

    sc = pcset.SetClass()  # A set-class for calculating names, etc.
    final_slices = []      # Holds the finalized slices to return
    first_measure_num_analyzed = -1     # We assume that the first measure is -1
    last_measure_num_analyzed = -1      # We assume that the last measure is -1
    current_measure_indices = [0 for i in range(len(parts))]  # The index of the next measure, for each part
    current_measure_num = -1      # The number of the next measure
    tempo = Decimal(60)    # We assume a tempo of 60 to begin
    time_signature = None  # The current time signature
    transpose = [0 for i in range(len(parts))]  # The amount by which to transpose, for each part

    # Adjust starting and ending measure numbers
    if starting_measure_num is None:
        starting_measure_num = -1
    if ending_measure_num is None:
        ending_measure_num = -1

    if len(parts) == 0:
        print("No parts were provided")

    else:
        # Determine the index of the first measure in each part. a is the part index,
        # and b is the index of the item inside the current part (which may or may not be a measure)
        for part_idx, part in enumerate(parts):
            found_first_measure = False
            for item_idx, item in enumerate(part):
                if type(item) == music21.stream.Measure:
                    if item.number >= starting_measure_num:
                        current_measure_num = item.number
                        first_measure_num_analyzed = current_measure_num
                        current_measure_indices[part_idx] = item_idx
                        found_first_measure = True
                # No need to continue after the first measure was found
                if found_first_measure:
                    break

        # We consider each measure separately. When we have finished the last measure,
        # the next_measure will reset to -1 and we will stop.
        while current_measure_num != -1:
            last_measure_num_analyzed = current_measure_num
            # The slices taken for this measure
            measure_slices = []

            # Consider each part separately for this measure
            for part_idx, part in enumerate(parts):
                # Tracks the number of slices taken for the current part in the current measure
                num_slices_taken = 0
                
                for item in part[current_measure_indices[part_idx]]:
                    last_item_was_voice = False
                    furthest_voice_slice = 0

                    # MusicXML doesn't handle transposition properly for 8va and 8vb clefs, so we need manual
                    # transposition. Record for the future.
                    if type(item) == music21.clef.Bass8vaClef or type(item) == music21.clef.Treble8vaClef:
                        transpose[part_idx] = 12
                    elif type(item) == music21.clef.Bass8vbClef or type(item) == music21.clef.Treble8vbClef:
                        transpose[part_idx] = -12
                    elif isinstance(item, music21.clef.Clef):
                        transpose[part_idx] = 0

                    # Record the current time signature
                    if type(item) == music21.meter.TimeSignature:
                        time_signature = item

                    # Update the tempo if we find a new one
                    # If a tempo map was provided, we check to see if the tempo changes in this measure
                    # If there was no tempo map, we check to see if an updated tempo exists in this measure
                    if tempo_map is not None:
                        if current_measure_num in tempo_map:
                            if type(tempo_map[current_measure_num]) == Fraction:
                                tempo = Decimal(tempo_map[current_measure_num].numerator) / Decimal(tempo_map[current_measure_num].denominator)
                            elif type(tempo_map[current_measure_num]) == Decimal:
                                tempo = tempo_map[current_measure_num]
                            else:
                                tempo = Decimal(tempo_map[current_measure_num])
                    elif type(item) == music21.tempo.MetronomeMark:
                        if item.number is not None:
                            tempo = Decimal(item.number)

                    # If we have found multiple voices in the same part in the same measure, we need to
                    # deal with each voice as a separate part.
                    if type(item) == music21.stream.Voice:
                        last_item_was_voice = True

                        # Track the start point for the voice
                        slice_start = num_slices_taken

                        for item2 in item:
                            # We can only take slices of notes, rests, or chords
                            if type(item2) == music21.note.Note or type(item2) == music21.note.Rest or type(
                                    item2) == music21.chord.Chord:
                                ql = item2.duration.quarterLength
                                if type(item2.duration.quarterLength) != Fraction:
                                    ql = Fraction(item2.duration.quarterLength)

                                # How many salami slices to take
                                num_slices = int(ql * n)
                                
                                # the pitches are considered as integers in p-space. 
                                # p_names_in_item hold pitch names which is often more convenient for humans.
                                pitches_in_item = []
                                p_names_in_item = []

                                # We use Morris's p-space. Obviously rests do not have pitches.
                                if type(item2) != music21.note.Rest:
                                    for p in item2.pitches:
                                        name = p.name
                                        octave = p.octave + (transpose[part_idx] // 12)
                                        pitches_in_item.append(p.midi - 60 + transpose[part_idx])
                                        p_names_in_item.append(name + str(octave))

                                # Perform salami slicing. num_slices is the number of slices we take for the current object.
                                for j in range(num_slices):
                                    if num_slices_taken >= len(measure_slices):
                                        measure_slices.append(
                                            SalamiSlice(tempo, Fraction(1, n.numerator), parts[part_idx][current_measure_indices[part_idx]].number,
                                                   len(parts)))
                                    measure_slices[num_slices_taken].add_pitches(pitches_in_item, p_names_in_item, part_idx)
                                    measure_slices[num_slices_taken].time_signature = time_signature
                                    measure_slices[num_slices_taken].start_position = Fraction(num_slices_taken,
                                                                                               n.numerator)
                                    num_slices_taken += 1

                        # Record the furthest slice reached in this voice if necessary
                        if furthest_voice_slice < num_slices_taken:
                            furthest_voice_slice = num_slices_taken

                        # Reset the slice counter to start on the next voice
                        num_slices_taken = slice_start

                    # If we just evaluated a voice and are done, we need to reset the slice counter
                    elif last_item_was_voice and num_slices_taken < furthest_voice_slice:
                        num_slices_taken = furthest_voice_slice

                    # We can only take slices of notes, rests, or chords. We cannot take slices of voices here;
                    # they are dealt with above.
                    if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(
                            item) == music21.chord.Chord:
                        ql = item.duration.quarterLength
                        if type(item.duration.quarterLength) != Fraction:
                            ql = Fraction(item.duration.quarterLength)

                        # How many salami slices to take
                        num_slices = int(ql * n)

                        # the pitches are considered as integers in p-space. 
                        # p_names_in_item holds pitch names which is often more convenient for humans.
                        pitches_in_item = []
                        p_names_in_item = []

                        # We use Morris's p-space. Obviously rests do not have pitches.
                        if type(item) != music21.note.Rest:
                            for p in item.pitches:
                                name = p.name
                                octave = p.octave + (transpose[part_idx] // 12)
                                pitches_in_item.append(p.midi - 60 + transpose[part_idx])
                                p_names_in_item.append(name + str(octave))

                        # Perform salami slicing. num_slices is the number of slices we take for the current object.
                        for j in range(num_slices):
                            if num_slices_taken >= len(measure_slices):
                                measure_slices.append(
                                    SalamiSlice(tempo, Fraction(1, n.numerator), parts[part_idx][current_measure_indices[part_idx]].number,
                                           len(parts)))
                            measure_slices[num_slices_taken].add_pitches(pitches_in_item, p_names_in_item, part_idx)
                            measure_slices[num_slices_taken].time_signature = time_signature
                            measure_slices[num_slices_taken].start_position = Fraction(num_slices_taken,
                                                                                       n.numerator)
                            num_slices_taken += 1

            # Clean up the slices from this measure
            clean_slices(measure_slices, True)
            for item in measure_slices:
                final_slices.append(item)

            # Find the next measure for each part
            for part_idx in range(len(parts)):  # a is the part index
                found_next = False
                current_measure_num = -1
                # We start at the item after the current measure
                for item_idx in range(current_measure_indices[part_idx] + 1, len(parts[part_idx])):
                    if type(parts[part_idx][item_idx]) == music21.stream.Measure:
                        current_measure_num = parts[part_idx][item_idx].number
                        current_measure_indices[part_idx] = item_idx
                        found_next = True
                    # No need to continue after the first measure was found
                    if found_next:
                        break

            # If we've analyzed the last measure, it's time to stop analyzing
            if current_measure_num > ending_measure_num > -1:
                current_measure_num = -1

    results = []
    global_bounds = get_piece_bounds(parts)

    # Make pctheory objects
    for sl in final_slices:
        sl.prepare_for_clean()

    clean_slices(final_slices, True, [section_divisions[i][0] for i in range(len(section_divisions))])
    for s in final_slices:
        s.run_calculations(sc)

    clean_slices(final_slices, False, [section_divisions[i][0] for i in range(len(section_divisions))])

    # Create sectional results
    for i in range(len(section_divisions)):
        section_slices = []
        start_time = 0
        for sl in final_slices:
            if sl.measure < section_divisions[i][0]:
                start_time += sl.duration
            elif sl.measure <= section_divisions[i][1]:
                section_slices.append(sl)
        bounds = global_bounds
        if use_local[i]:
            bounds = get_bounds(section_slices)
        set_slice_bounds(section_slices, bounds)
        for s in section_slices:
            s.run_calculations_burt()
        results.append(Results(section_slices, section_divisions[i][0], section_divisions[i][1],
                               len(parts), start_time))

    # Create overall results
    clean_slices(final_slices)
    bounds = global_bounds
    if len(use_local) == 1:
        if use_local[0]:
            bounds = get_bounds(final_slices)
    set_slice_bounds(final_slices, bounds)
    for f_slice in final_slices:
        f_slice.run_calculations_burt()
    results.insert(0, Results(final_slices, first_measure_num_analyzed, last_measure_num_analyzed, len(parts)))
    return results


def read_analysis_from_file(path):
    """
    Reads analysis data from a JSON file
    :param path: The file path
    :return: A list of Results objects
    """
    data = None
    results = []
    with open(path, "r") as file_in:
        data = json.load(file_in)
    for item in data:
        slices = []
        for dslice in item["slices"]:
            cslice = SalamiSlice()
            cslice._chord_spacing_contour = dslice["chord_spacing_contour"]
            cslice._chord_spacing_index = float(dslice["chord_spacing_index"])
            cslice._core = bool(dslice["core"])
            cslice._derived_core = bool(dslice["derived_core"])
            cslice._derived_core_associations = dslice["derived_core_associations"]
            cslice._duration = Decimal(dslice["duration"])
            if dslice["ioi_in_seconds"] != "None":
                cslice._ioi_in_seconds = float(dslice["ioi_in_seconds"])
            else:
                cslice._ioi_in_seconds = None
            cslice._ipseg = dslice["ipseg"]
            cslice._measure = dslice["measure"]
            cslice._pset_cardinality = dslice["pset_cardinality"]
            cslice._pitch_count_with_duplicates = dslice["pitch_count_with_duplicates"]
            cslice._pcset_cardinality = dslice["pcset_cardinality"]
            cslice._pcseg = [pitch.PitchClass(pc) for pc in dslice["pcseg"]]
            cslice._pcset = set(cslice.pcseg)
            cslice._pcsegs = [[pitch.PitchClass(pc) for pc in dslice["pcsegs"][v]]
                              for v in range(len(dslice["pcsegs"]))]
            cslice._pcset = set(cslice.pcseg)
            cslice._pcsets = [set(cslice.pcsegs[v]) for v in range(len(cslice.pcsegs))]
            cslice._pitchseg = dslice["pitchseg"]
            cslice._pitchsegs = dslice["pitchsegs"]
            cslice._pitch_name_list = dslice["pitch_name_list"]
            cslice._pitch_name_lists = dslice["pitch_name_lists"]
            cslice._pseg = [pitch.Pitch(p) for p in dslice["pseg"]]
            cslice._psegs = [[pitch.Pitch(p) for p in dslice["psegs"][v]] for v in range(len(dslice["psegs"]))]
            cslice._pset = set(cslice.pseg)
            cslice._psets = [set(cslice.psegs[v]) for v in range(len(cslice.psegs))]
            cslice._quarter_duration = Fraction(dslice["quarter_duration"][0], dslice["quarter_duration"][1])
            cslice._sc_name = dslice["sc_name"]
            cslice._sc_name_carter = dslice["sc_name_carter"]
            cslice._ins = dslice["ins"]
            cslice._lns = dslice["lns"]
            cslice._lower_bound = dslice["lower_bound"]
            cslice._median_trajectory = dslice["median_trajectory"]
            cslice._ns = dslice["ns"]
            cslice._ps = dslice["ps"]
            cslice._start_position = Fraction(dslice["start_position"][0], dslice["start_position"][1])
            cslice._time_signature = music21.meter.TimeSignature(dslice["time_signature"])
            cslice._uns = dslice["uns"]
            cslice._upper_bound = dslice["upper_bound"]
            slices.append(cslice)
        result = Results(slices, item["measure_num_first"], item["measure_num_last"], len(item["pitch_highest_voices"]))
        result._max_pitch_count_with_duplicates = item["max_pitch_count_with_duplicates"]
        result._chord_spacing_contour_duration = {}
        result._chord_spacing_contour_frequency = item["chord_spacing_contour_frequency"]
        result._duration = Decimal(item["duration"])
        result._duration_avg = Decimal(item["duration_avg"])
        result._ins_avg = item["ins_avg"]
        result._ins_max = item["ins_max"]
        result._ins_min = item["ins_min"]
        result._ioi_avg_in_seconds = Decimal(item["ioi_avg_in_seconds"])
        result._lns_avg = item["lns_avg"]
        result._lns_max = item["lns_max"]
        result._lns_min = item["lns_min"]
        result._lower_bound = item["lower_bound"]
        result._lps_card = item["lps_card"]
        result._median_trajectory_avg = item["median_trajectory_avg"]
        result._median_trajectory_max = item["median_trajectory_max"]
        result._median_trajectory_min = item["median_trajectory_min"]
        result._num_measures = item["num_measures"]
        result._num_voices = len(item["pitch_highest_voices"])
        result._pitch_highest = item["pitch_highest"]
        result._pitch_highest_voices = item["pitch_highest_voices"]
        result._pitch_lowest = item["pitch_lowest"]
        result._pitch_lowest_voices = item["pitch_lowest_voices"]
        result._pset_card_avg = item["pset_card_avg"]
        result._chord_spacing_index_avg = item["chord_spacing_index_avg"]
        result._pset_duration = {}
        result._pset_frequency = item["pset_frequency"]
        result._psc_duration = {}
        result._pcsc_frequency = item["pcsc_frequency"]
        result._psc_frequency = item["psc_frequency"]
        result._ps_avg = item["ps_avg"]
        result._ps_max = item["ps_max"]
        result._ps_min = item["ps_min"]
        result._quarter_duration = Fraction(item["quarter_duration"][0], item["quarter_duration"][1])
        result._start_time = Decimal(item["start_time"])
        result._uns_avg = item["uns_avg"]
        result._uns_max = item["uns_max"]
        result._uns_min = item["uns_min"]
        result._upper_bound = item["upper_bound"]
        result._pc_duration = {}
        result._pc_frequency = {}
        result._pitch_duration = {}
        result._pitch_frequency = {}
        for key, val in item["chord_spacing_contour_duration"].items():
            result.chord_spacing_contour_duration[key] = Decimal(val)
        for key, val in item["pc_duration"].items():
            result.pc_duration[int(key)] = Decimal(val)
        for v in range(len(item["pc_duration_voices"])):
            result.pc_duration_voices.append({})
            for key, val in item["pc_duration_voices"][v].items():
                result.pc_duration_voices[v][int(key)] = Decimal(val)
        for key, val in item["pitch_duration"].items():
            result.pitch_duration[int(key)] = Decimal(val)
        for v in range(len(item["pitch_duration_voices"])):
            result.pitch_duration_voices.append({})
            for key, val in item["pitch_duration_voices"][v].items():
                result.pitch_duration_voices[v][int(key)] = Decimal(val)
        for key, val in item["pc_frequency"].items():
            result.pc_frequency[int(key)] = val
        for v in range(len(item["pc_frequency_voices"])):
            result.pc_frequency_voices.append({})
            for key, val in item["pc_frequency_voices"][v].items():
                result.pc_frequency_voices[v][int(key)] = val
        for key, val in item["pitch_frequency"].items():
            result.pitch_frequency[int(key)] = val
        for v in range(len(item["pitch_frequency_voices"])):
            result.pitch_frequency_voices.append({})
            for key, val in item["pitch_frequency_voices"][v].items():
                result.pitch_frequency_voices[v][int(key)] = val
        for key, val in item["pset_duration"].items():
            result.pset_duration[key] = Decimal(val)
        for key, val in item["pcsc_duration"].items():
            result.pcsc_duration[key] = Decimal(val)
        for key, val in item["psc_duration"].items():
            result.psc_duration[key] = Decimal(val)
        results.append(result)
    return results


def write_analysis_to_file(results, path):
    """
    Writes an analysis to a JSON file
    :param results: A Results object
    :param path: A path for a JSON file
    :return:
    """
    data = []
    for i in range(len(results)):
        data.append({})
        data[i]["max_pitch_count_with_duplicates"] = results[i].max_pitch_count_with_duplicates
        data[i]["chord_spacing_contour_duration"] = {}
        data[i]["chord_spacing_contour_frequency"] = results[i].chord_spacing_contour_frequency
        data[i]["duration"] = str(results[i].duration)
        data[i]["duration_avg"] = str(results[i].duration_avg)
        data[i]["ins_avg"] = results[i].ins_avg
        data[i]["ins_max"] = results[i].ins_max
        data[i]["ins_min"] = results[i].ins_min
        data[i]["ioi_avg_in_seconds"] = str(results[i].ioi_avg_in_seconds)
        data[i]["lns_avg"] = results[i].lns_avg
        data[i]["lns_max"] = results[i].lns_max
        data[i]["lns_min"] = results[i].lns_min
        data[i]["lower_bound"] = results[i].lower_bound
        data[i]["lps_card"] = results[i].lps_card
        data[i]["measure_num_first"] = results[i].measure_num_first
        data[i]["measure_num_last"] = results[i].measure_num_last
        data[i]["median_trajectory_avg"] = results[i].median_trajectory_avg
        data[i]["median_trajectory_max"] = results[i].median_trajectory_max
        data[i]["median_trajectory_min"] = results[i].median_trajectory_min
        data[i]["num_measures"] = results[i].num_measures
        data[i]["pitch_highest"] = results[i].pitch_highest
        data[i]["pitch_highest_voices"] = results[i].pitch_highest_voices
        data[i]["pitch_lowest"] = results[i].pitch_lowest
        data[i]["pitch_lowest_voices"] = results[i].pitch_lowest_voices
        data[i]["pset_card_avg"] = results[i].pset_card_avg
        data[i]["chord_spacing_index_avg"] = results[i].chord_spacing_index_avg
        data[i]["pset_duration"] = {}
        data[i]["pset_frequency"] = results[i].pset_frequency
        data[i]["pcsc_duration"] = {}
        data[i]["psc_duration"] = {}
        data[i]["pcsc_frequency"] = results[i].pcsc_frequency
        data[i]["psc_frequency"] = results[i].psc_frequency
        data[i]["ps_avg"] = results[i].ps_avg
        data[i]["ps_max"] = results[i].ps_max
        data[i]["ps_min"] = results[i].ps_min
        data[i]["quarter_duration"] = [results[i].quarter_duration.numerator, results[i].quarter_duration.denominator]
        data[i]["start_time"] = str(results[i].start_time)
        data[i]["uns_avg"] = results[i].uns_avg
        data[i]["uns_max"] = results[i].uns_max
        data[i]["uns_min"] = results[i].uns_min
        data[i]["upper_bound"] = results[i].upper_bound
        data[i]["pc_duration"] = {}
        data[i]["pc_duration_voices"] = []
        data[i]["pc_frequency"] = results[i].pc_frequency
        data[i]["pc_frequency_voices"] = results[i].pc_frequency_voices
        data[i]["pitch_duration"] = {}
        data[i]["pitch_duration_voices"] = []
        data[i]["pitch_frequency"] = results[i].pitch_frequency
        data[i]["pitch_frequency_voices"] = results[i].pitch_frequency_voices
        data[i]["slices"] = []
        for key, val in results[i].chord_spacing_contour_duration.items():
            data[i]["chord_spacing_contour_duration"][key] = str(val)
        for key, val in results[i].pc_duration.items():
            data[i]["pc_duration"][key] = str(val)
        for v in range(len(results[i].pc_duration_voices)):
            data[i]["pc_duration_voices"].append({})
            for key, val in results[i].pc_duration_voices[v].items():
                data[i]["pc_duration_voices"][len(data[i]["pc_duration_voices"]) - 1][key] = str(val)
        for key, val in results[i].pitch_duration.items():
            data[i]["pitch_duration"][key] = str(val)
        for v in range(len(results[i].pitch_duration_voices)):
            data[i]["pitch_duration_voices"].append({})
            for key, val in results[i].pitch_duration_voices[v].items():
                data[i]["pitch_duration_voices"][len(data[i]["pitch_duration_voices"]) - 1][key] = str(val)
        for key, val in results[i].pset_duration.items():
            data[i]["pset_duration"][key] = str(val)
        for key, val in results[i].pcsc_duration.items():
            data[i]["pcsc_duration"][key] = str(val)
        for key, val in results[i].psc_duration.items():
            data[i]["psc_duration"][key] = str(val)
        for rslice in results[i].slices:
            cslice = {}
            cslice["chord_spacing_contour"] = rslice.chord_spacing_contour
            cslice["chord_spacing_index"] = rslice.chord_spacing_index
            cslice["core"] = int(rslice.core)
            cslice["derived_core"] = int(rslice.derived_core)
            cslice["derived_core_associations"] = rslice.derived_core_associations
            cslice["duration"] = str(rslice.duration)
            cslice["ipseg"] = rslice.ipseg
            cslice["ioi_in_seconds"] = str(rslice.ioi_in_seconds)
            cslice["measure"] = rslice.measure
            cslice["pset_cardinality"] = rslice.pset_cardinality
            cslice["pitch_count_with_duplicates"] = rslice.pitch_count_with_duplicates
            cslice["pcset_cardinality"] = rslice.pcset_cardinality
            cslice["pcseg"] = [pc.pc for pc in rslice.pcseg]
            cslice["pcsegs"] = [[pc.pc for pc in rslice.pcsegs[v]] for v in range(len(rslice.pcsegs))]
            cslice["pitchseg"] = [p for p in rslice.pitchseg]
            cslice["pitchsegs"] = [[p for p in rslice.pitchsegs[v]] for v in range(len(rslice.pitchsegs))]
            cslice["pitch_name_list"] = [pname for pname in rslice.pitch_name_list]
            cslice["pitch_name_lists"] = [[pname for pname in rslice.pitch_name_lists[v]] for v in range(len(rslice.pitch_name_lists))]
            cslice["pseg"] = [p.p for p in rslice.pseg]
            cslice["psegs"] = [[p.p for p in rslice.psegs[v]] for v in range(len(rslice.psegs))]
            cslice["quarter_duration"] = [rslice.quarter_duration.numerator, rslice.quarter_duration.denominator]
            cslice["sc_name"] = rslice.sc_name
            cslice["sc_name_carter"] = rslice.sc_name_carter
            cslice["ins"] = rslice.ins
            cslice["lns"] = rslice.lns
            cslice["lower_bound"] = rslice.lower_bound
            cslice["median_trajectory"] = rslice.median_trajectory
            cslice["ns"] = rslice.ns
            cslice["ps"] = rslice.ps
            cslice["start_position"] = [rslice.start_position.numerator, rslice.start_position.denominator]
            cslice["time_signature"] = rslice.time_signature.ratioString
            cslice["uns"] = rslice.uns
            cslice["upper_bound"] = rslice.upper_bound
            data[i]["slices"].append(cslice)

    with open(path, "w") as out:
        out.write(json.dumps(data))


def write_general_report(section_name, file, file_command, results, lowest_pitch, highest_pitch):
    """
    Writes a general report to Excel
    :param section_name: The name of the section being reported
    :param file: The file path
    :param file_command: The command ("w" or "a")
    :param results: A Results object
    :param lowest_pitch: The lowest pitch analyzed
    :param highest_pitch: The highest pitch analyzed
    :return: None
    """
    report = {
        "section": [],
        "start_time": [], 
        "duration": [], 
        "duration_avg": [], 
        "ioi_avg": [], 
        "pset_card_avg": [],
        "pitch_count_with_duplicate_avg": [],
        "csi_avg": [], 
        "lps": [], 
        "p_u": [], 
        "p_l": [], 
        "ps_avg": [], 
        "ps_min": [], 
        "ps_max": [],
        "uns_avg": [], 
        "uns_min": [], 
        "uns_max": [], 
        "lns_avg": [], 
        "lns_min": [], 
        "lns_max": [], 
        "ins_avg": [], 
        "ins_min": [], 
        "ins_max": [], 
        "mt_avg": [], 
        "mt_min": [],
        "mt_max": []
    }
    for i in range(0, 12):
        report[f"pc{i}_dur"] = []
    for i in range(0, 12):
        report[f"pc{i}_freq"] = []
    for i in range(lowest_pitch, highest_pitch + 1):
        report[f"p{i}_dur"] = []
    for i in range(lowest_pitch, highest_pitch + 1):
        report[f"p{i}_freq"] = []

    report["section"].append(section_name)
    report["start_time"].append(float(results.start_time))
    report["duration"].append(float(results.duration))
    report["duration_avg"].append(float(results.duration_avg))
    report["ioi_avg"].append(float(results.ioi_avg_in_seconds))
    report["pset_card_avg"].append(results.pset_card_avg)
    report["pitch_count_with_duplicate_avg"].append(results.pitch_count_with_duplicates_avg)
    report["csi_avg"].append(results.chord_spacing_index_avg)
    report["lps"].append(results.lps_card)
    report["p_u"].append(results.pitch_highest)
    report["p_l"].append(results.pitch_lowest)
    report["ps_avg"].append(results.ps_avg)
    report["ps_min"].append(results.ps_min)
    report["ps_max"].append(results.ps_max)
    report["uns_avg"].append(results.uns_avg)
    report["uns_min"].append(results.uns_min)
    report["uns_max"].append(results.uns_max)
    report["lns_avg"].append(results.lns_avg)
    report["lns_min"].append(results.lns_min)
    report["lns_max"].append(results.lns_max)
    report["ins_avg"].append(results.ins_avg)
    report["ins_min"].append(results.ins_min)
    report["ins_max"].append(results.ins_max)
    report["mt_avg"].append(results.median_trajectory_avg)
    report["mt_min"].append(results.median_trajectory_min)
    report["mt_max"].append(results.median_trajectory_max)
    
    for i in range(0, 12):
        if i in results.pc_duration.keys():
            report[f"pc{i}_dur"].append(float(results.pc_duration[i]))
        else:
            report[f"pc{i}_dur"].append(0)
    for i in range(0, 12):
        if i in results.pc_frequency.keys():
            report[f"pc{i}_freq"].append(results.pc_frequency[i])
        else:
            report[f"pc{i}_freq"].append(0)
    for i in range(lowest_pitch, highest_pitch + 1):
        if i in results.pitch_duration.keys():
            report[f"p{i}_dur"].append(float(results.pitch_duration[i]))
        else:
            report[f"p{i}_dur"].append(0)
    for i in range(lowest_pitch, highest_pitch + 1):
        if i in results.pitch_frequency.keys():
            report[f"p{i}_freq"].append(results.pitch_frequency[i])
        else:
            report[f"p{i}_freq"].append(0)

    for v in range(results.num_voices):
        report["section"].append(f"{section_name} (Voice {v})")
        report["start_time"].append(None)
        report["duration"].append(None)
        report["duration_avg"].append(None)
        report["ioi_avg"].append(None)
        report["pset_card_avg"].append(None)
        report["pitch_count_with_duplicate_avg"].append(None)
        report["csi_avg"].append(None)
        report["lps"].append(results.pitch_highest_voices[v] - results.pitch_lowest_voices[v] + 1)
        report["p_u"].append(results.pitch_highest_voices[v])
        report["p_l"].append(results.pitch_lowest_voices[v])
        report["ps_avg"].append(None)
        report["ps_min"].append(None)
        report["ps_max"].append(None)
        report["uns_avg"].append(None)
        report["uns_min"].append(None)
        report["uns_max"].append(None)
        report["lns_avg"].append(None)
        report["lns_min"].append(None)
        report["lns_max"].append(None)
        report["ins_avg"].append(None)
        report["ins_min"].append(None)
        report["ins_max"].append(None)
        report["mt_avg"].append(None)
        report["mt_min"].append(None)
        report["mt_max"].append(None)
        for i in range(0, 12):
            if i in results.pc_duration_voices[v].keys():
                report[f"pc{i}_dur"].append(float(results.pc_duration_voices[v][i]))
            else:
                report[f"pc{i}_dur"].append(0)
        for i in range(0, 12):
            if i in results.pc_frequency_voices[v].keys():
                report[f"pc{i}_freq"].append(results.pc_frequency_voices[v][i])
            else:
                report[f"pc{i}_freq"].append(0)
        for i in range(lowest_pitch, highest_pitch + 1):
            if i in results.pitch_duration_voices[v].keys():
                report[f"p{i}_dur"].append(float(results.pitch_duration_voices[v][i]))
            else:
                report[f"p{i}_dur"].append(0)
        for i in range(lowest_pitch, highest_pitch + 1):
            if i in results.pitch_frequency_voices[v].keys():
                report[f"p{i}_freq"].append(results.pitch_frequency_voices[v][i])
            else:
                report[f"p{i}_freq"].append(0)
        
    output_df = pd.DataFrame(report)
    if file_command == 'a':
        contents = pd.read_excel(file)
        new_df = pd.concat([contents, output_df])
        new_df.to_excel(file, index=False)
    else:
        output_df.to_excel(file, index=False)


def write_statistics(file, headings, dictionaries):
    """
    Writes a dictionary to Excel
    :param file: A file name
    :param headings: A headings row for the file
    :param dictionaries: Dictionaries with common keys to write to file
    :return: None
    """
    stat_list = []
    for key, value in dictionaries[0].items():
        stat_list.append([key, value])
    for i in range(1, len(dictionaries)):
        for j in range(len(stat_list)):
            stat_list[j].append(dictionaries[i][stat_list[j][0]])
    stat_list = sorted(stat_list, key=lambda x: x[0])
    stat_list = sorted(stat_list, key=lambda x: len(x[0]))

    output = [[] for heading in headings]
    for line in stat_list:
        if len(line) > 0:
            output[0].append(f"{line[0]}")
        for i in range(1, len(line)):
            if type(line[i]) == Decimal:
                output[i].append(float(line[i]))
            else:
                output[i].append(line[i])
    output_df = pd.DataFrame({headings[i]: output[i] for i in range(len(output))})
    output_df.to_excel(file, index=False)


def write_report(file, results):
    """
    Writes a report to Excel
    :param file: A file name (and path if necessary)
    :param results: A Results object
    """
    report = {
        "measure_no": [],
        "start_time": [],
        "duration": [],
        "ioi": [],
        "quarter_length": [],
        "chord_cardinality": [],
        "ps": [],
        "match": [],
        "ns": [],
        "uns": [],
        "ins": [],
        "lns": [],
        "mt": [],
        "morris_name": [],
        "carter_name": [],
        "core": [],
        "derived_core": [],
        "dc_associations": [],
        "pcset": [],
        "pset": [],
        "psc": [],
        "chord_spacing_contour": [],
        "chord_spacing_index": []
    }
    for i in range(results.max_pitch_count_with_duplicates):
        report[f"pitch_{i + 1}"] = []
    for i in range(results.ps_max):
        report[f"pn_{i + 1}"] = []

    if len(results.slices) > 0:
        # Track the onset position in seconds
        position = 0

        # Output each slice
        for item in results.slices:
            report["measure_no"].append(item.measure)
            report["start_time"].append(float(position))
            report["duration"].append(float(item.duration))
            if item.ioi_in_seconds is None:
                report["ioi"].append(item.ioi_in_seconds)
            else:
                report["ioi"].append(float(item.ioi_in_seconds))
            report["quarter_length"].append(item.quarter_duration)
            report["chord_cardinality"].append(item.pitch_count_with_duplicates)
            report["ps"].append(item.ps)
            report["match"].append(item.pitch_count_with_duplicates == item.ps)
            report["ns"].append(item.ns)
            report["uns"].append(item.uns)
            report["ins"].append(item.ins)
            report["lns"].append(item.lns)
            report["mt"].append(item.median_trajectory)
            report["morris_name"].append(item.sc_name)
            report["carter_name"].append(item.sc_name_carter)
            report["core"].append(item.core)
            report["derived_core"].append(item.derived_core)
            report["dc_associations"].append(item.derived_core_associations)
            report["pcset"].append(item.get_pcset_string())
            report["pset"].append(item.get_pset_string())
            report["psc"].append(item.get_ipseg_string())
            report["chord_spacing_contour"].append(item.get_chord_spacing_contour_string())
            report["chord_spacing_index"].append(item.chord_spacing_index)
            for i in range(results.max_pitch_count_with_duplicates):
                if i < len(item.pitchseg):
                    report[f"pitch_{i + 1}"].append(item.pitch_name_list[i])
                else:
                    report[f"pitch_{i + 1}"].append(None)
            for i in range(results.ps_max):
                if i < len(item.pseg):
                    report[f"pn_{i + 1}"].append(item.pseg[i])
                else:
                    report[f"pn_{i + 1}"].append(None)
            position += item.duration

    report_df = pd.DataFrame(report)
    report_df.to_excel(file, index=False)
