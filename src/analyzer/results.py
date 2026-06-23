"""
File: results.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains the Results class for analysis results.
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

import numpy as np
from decimal import Decimal


class Results:
    def __init__(self, slices, measure_num_first, measure_num_last, voices, start_time=0):
        """
        Creates a Results object
        :param slices: A list of slices
        :param measure_num_first: The first measure number analyzed
        :param measure_num_last: The last measure number analyzed
        """
        self._max_pitch_count_with_duplicates = 0  # The maximum number of pitches in a chord (may be greater than PS)
        self._chord_spacing_contour_duration = None
        self._chord_spacing_contour_frequency = None
        self._duration = 0
        self._duration_avg = 0  # The slice duration average. Note that this is not IOI; it includes slices that are rests.
        self._ins_avg = 0       # The INS average
        self._ins_max = 0
        self._ins_min = 0
        self._ioi_avg_in_seconds = 0       # The IOI (interonset interval) average
        self._lns_avg = 0       # The LNS average
        self._lns_max = 0
        self._lns_min = 0
        self._lower_bound = 0
        self._lps_card = 0      # The cardinality of LPS
        self._measure_num_first = measure_num_first
        self._measure_num_last = measure_num_last
        self._median_trajectory_avg = 0   # The MT average
        self._median_trajectory_max = 0
        self._median_trajectory_min = 0
        self._num_measures = measure_num_last - measure_num_first + 1
        self._num_voices = voices
        self._pc_duration = None
        self._pc_duration_voices = None
        self._pc_frequency = None
        self._pc_frequency_voices = None
        self._pitch_duration = None
        self._pitch_duration_voices = None
        self._pitch_frequency = None
        self._pitch_frequency_voices = None
        self._pitch_highest = -np.inf
        self._pitch_highest_voices = None
        self._pitch_lowest = np.inf
        self._pitch_lowest_voices = None
        self._pset_card_avg = 0
        self._pitch_count_with_duplicates_avg = 0
        self._pset_duration = None
        self._pset_frequency = None
        self._chord_spacing_index_avg = 0
        self._pcsc_duration = None
        self._pcsc_frequency = None
        self._psc_duration = None
        self._psc_frequency = None
        self._ps_avg = 0  # The PS average
        self._ps_max = 0
        self._ps_min = 0
        self._quarter_duration = 0
        self._slices = slices
        self._start_time = Decimal(start_time)
        self._uns_avg = 0  # The UNS average
        self._uns_max = 0
        self._uns_min = 0
        self._upper_bound = 0

        # Calculate values
        self._calculate_values()

    @property
    def max_pitch_count_with_duplicates(self):
        """
        The maximum cardinality of the analyzed SalamiSlices
        :return: The maximum cardinality of the analyzed SalamiSlices
        """
        return self._max_pitch_count_with_duplicates

    @property
    def chord_spacing_contour_duration(self):
        """
        A dictionary in which the chord spacing contours that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._chord_spacing_contour_duration

    @property
    def chord_spacing_contour_frequency(self):
        """
        A dictionary in which the chord spacing contours that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._chord_spacing_contour_frequency

    @property
    def duration(self):
        """
        The total duration of all the SalamiSlices in seconds
        :return: The total duration of all the SalamiSlices in seconds
        """
        return self._duration

    @property
    def duration_avg(self):
        """
        The average duration of all the SalamiSlices in seconds
        :return: The average duration of all the SalamiSlices in seconds
        """
        return self._duration_avg

    @property
    def ins_avg(self):
        """
        The average internal negative space (INS)
        :return: The average internal negative space (INS)
        """
        return self._ins_avg

    @property
    def ins_max(self):
        """
        The maximum internal negative space (INS)
        :return: The maximum internal negative space (INS)
        """
        return self._ins_max

    @property
    def ins_min(self):
        """
        The minimum internal negative space (INS)
        :return: The minimum internal negative space (INS)
        """
        return self._ins_min
    
    @property
    def ioi_avg_in_seconds(self):
        """
        The average interonset interval (IOI) in *seconds*
        :return: The average interonset interval (IOI) in *seconds*
        """
        return self._ioi_avg_in_seconds
    
    @property
    def lns_avg(self):
        """
        The average lower negative space (LNS)
        :return: The average LNS
        """
        return self._lns_avg

    @property
    def lns_max(self):
        """
        The maximum lower negative space (LNS)
        :return: The maximum lower negative space (LNS)
        """
        return self._lns_max

    @property
    def lns_min(self):
        """
        The minimum lower negative space (LNS)
        :return: The minimum lower negative space (LNS)
        """
        return self._lns_min

    @property
    def lower_bound(self):
        """
        The lower bound of the SalamiSlices
        :return: The lower bound of the SalamiSlices
        """
        return self._lower_bound

    @property
    def lps_card(self):
        """
        The cardinality of local pitch space (LPS)
        :return: The cardinality of local pitch space (LPS)
        """
        return self._lps_card

    @property
    def measure_num_first(self):
        """
        The first measure number analyzed
        :return: The first measure number analyzed
        """
        return self._measure_num_first

    @property
    def measure_num_last(self):
        """
        The last measure number analyzed
        :return: The last measure number analyzed
        """
        return self._measure_num_last

    @property
    def median_trajectory_avg(self):
        """
        The average median trajectory
        :return: The average median trajectory
        """
        return self._median_trajectory_avg

    @property
    def median_trajectory_max(self):
        """
        The maximum median trajectory
        :return: The maximum median trajectory
        """
        return self._median_trajectory_max

    @property
    def median_trajectory_min(self):
        """
        The minimum median trajectory
        :return: The minimum median trajectory
        """
        return self._median_trajectory_min

    @property
    def num_measures(self):
        """
        The number of measures analyzed
        :return: The number of measures analyzed
        """
        return self._num_measures

    @property
    def num_voices(self):
        """
        The number of voices
        :return: The number of voices
        """
        return self._num_voices

    @property
    def pc_duration(self):
        """
        A dictionary in which the pcs that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pc_duration

    @property
    def pc_frequency(self):
        """
        A dictionary in which the pcs that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pc_frequency

    @property
    def pitch_duration(self):
        """
        A dictionary in which the pitches that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pitch_duration

    @property
    def pitch_frequency(self):
        """
        A dictionary in which the pitches that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pitch_frequency

    @property
    def pc_duration_voices(self):
        """
        A dictionary in which the pcs that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pc_duration_voices

    @property
    def pc_frequency_voices(self):
        """
        A dictionary in which the pcs that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pc_frequency_voices

    @property
    def pitch_duration_voices(self):
        """
        A dictionary in which the pitches that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pitch_duration_voices

    @property
    def pitch_frequency_voices(self):
        """
        A dictionary in which the pitches that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pitch_frequency_voices

    @property
    def pitch_highest(self):
        """
        The highest pitch in the analyzed measures
        :return: The highest pitch in the analyzed measures
        """
        return self._pitch_highest

    @property
    def pitch_lowest(self):
        """
        The lowest pitch in the analyzed measures
        :return: The lowest pitch in the analyzed measures
        """
        return self._pitch_lowest

    @property
    def pitch_highest_voices(self):
        """
        The highest pitch in the analyzed measures
        :return: The highest pitch in the analyzed measures
        """
        return self._pitch_highest_voices

    @property
    def pitch_lowest_voices(self):
        """
        The lowest pitch in the analyzed measures
        :return: The lowest pitch in the analyzed measures
        """
        return self._pitch_lowest_voices

    @property
    def pset_card_avg(self):
        """
        The average pset cardinality of the analyzed measures (by duration)
        :return: The average pset cardinality
        """
        return self._pset_card_avg
    
    @property
    def pitch_count_with_duplicates_avg(self):
        """
        The average pitch count with duplicates of the analyzed measures (by duration).
        Note that this is very similar to `pset_card_avg`, except that this metric includes
        duplicate pitches in its calculation.
        :return: The average pitch count with duplicates
        """
        return self._pitch_count_with_duplicates_avg

    @property
    def chord_spacing_index_avg(self):
        """
        The average pset spacing index of the analyzed measures (by duration)
        :return: The average pset spacing index
        """
        return self._chord_spacing_index_avg

    @property
    def pset_duration(self):
        """
        A dictionary in which the psets that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pset_duration

    @property
    def pset_frequency(self):
        """
        A dictionary in which the psets that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pset_frequency

    @property
    def pcsc_duration(self):
        """
        A dictionary in which the pscs that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._pcsc_duration

    @property
    def pcsc_frequency(self):
        """
        A dictionary in which the pscs that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._pcsc_frequency

    @property
    def psc_duration(self):
        """
        A dictionary in which the pscs that occur in the analyzed measures are the keys,
        and their cumulative durations in the analyzed measures (in seconds) are the values
        :return: A dictionary
        """
        return self._psc_duration

    @property
    def psc_frequency(self):
        """
        A dictionary in which the pscs that occur in the analyzed measures are the keys,
        and the number of nonconsecutive occurrences in the analyzed measures are the values
        :return: A dictionary
        """
        return self._psc_frequency

    @property
    def ps_avg(self):
        """
        The average positive space (PS)
        :return: The average positive space (PS)
        """
        return self._ps_avg

    @property
    def ps_max(self):
        """
        The maximum positive space (PS)
        :return: The maximum positive space (PS)
        """
        return self._ps_max

    @property
    def ps_min(self):
        """
        The minimum positive space (PS)
        :return: The minimum positive space (PS)
        """
        return self._ps_min

    @property
    def quarter_duration(self):
        """
        The duration of the SalamiSlices, in quarter notes
        :return: The duration of the SalamiSlices, in quarter notes
        """
        return self._quarter_duration

    @property
    def slices(self):
        """
        The SalamiSlices
        :return: The SalamiSlices
        """
        return self._slices

    @property
    def start_time(self):
        """
        The start time of the Results object
        :return: The start time
        """
        return self._start_time

    @property
    def uns_avg(self):
        """
        The average upper negative space (UNS)
        :return: The average upper negative space (UNS)
        """
        return self._uns_avg

    @property
    def uns_max(self):
        """
        The maximum upper negative space (UNS)
        :return: The maximum upper negative space (UNS)
        """
        return self._uns_max

    @property
    def uns_min(self):
        """
        The minimum upper negative space (UNS)
        :return: The minimum upper negative space (UNS)
        """
        return self._uns_min

    @property
    def upper_bound(self):
        """
        The upper bound of the SalamiSlices
        :return: The upper bound of the SalamiSlices
        """
        return self._upper_bound

    def _calculate_values(self):
        """
        Calculates values for the results object
        :return: None
        """
        valid_durations_for_spacing = 0
        if len(self._slices) > 0:
            # Establish basic information
            self._lower_bound = 0
            self._lps_card = 0
            self._upper_bound = 0
            self._ins_min = 0
            self._lns_min = 0
            self._median_trajectory_max = 0
            self._median_trajectory_min = 0
            
            i = 0
            while i < len(self._slices):
                if self._slices[i].lower_bound is not None and self._slices[i].upper_bound is not None:
                    self._lower_bound = self._slices[0].lower_bound
                    self._lps_card = self._slices[0].upper_bound - self._slices[0].lower_bound + 1
                    self._upper_bound = self._slices[0].upper_bound
                    self._ins_min = self._lps_card
                    self._lns_min = self._lps_card
                    self._median_trajectory_max = self._lower_bound
                    self._median_trajectory_min = self._upper_bound
                    break
                i += 1
            
            self._chord_spacing_contour_duration = {}  # The duration of each cseg
            self._chord_spacing_contour_frequency = {}  # The number of occurrences of each cseg
            self._pc_duration = {}  # The total number of seconds that this pitch-class is active
            self._pc_duration_voices = [{} for v in range(self._num_voices)]
            self._pc_frequency = {}  # The total number of distinct (nonadjacent) occurrences of this pitch-class
            self._pc_frequency_voices = [{} for v in range(self._num_voices)]
            self._pitch_duration = {}  # The total number of seconds that this pitch is active
            self._pitch_duration_voices = [{} for v in range(self._num_voices)]
            self._pitch_frequency = {}  # The total number of distinct (nonadjacent) occurrences of this pitch
            self._pitch_frequency_voices = [{} for v in range(self._num_voices)]
            self._pitch_highest_voices = [-np.inf for v in range(self._num_voices)]
            self._pitch_lowest_voices = [np.inf for v in range(self._num_voices)]
            self._pset_duration = {}
            self._pset_frequency = {}
            self._pcsc_duration = {}
            self._pcsc_frequency = {}
            self._psc_duration = {}
            self._psc_frequency = {}
            self._ps_min = self._lps_card
            self._uns_min = self._lps_card
            valid_durations_for_spacing = 0
            
            ioi_count = 0

            # Sum values for averages and establish maxes and mins
            for s in self._slices:
                self._duration += s.duration
                self._quarter_duration += s.quarter_duration
                if s.pset_cardinality > 0:
                    self._pset_card_avg += len(s.pset) * s.duration
                    self._pitch_count_with_duplicates_avg += s.pitch_count_with_duplicates * s.duration
                    self._ioi_avg_in_seconds += s.ioi_in_seconds
                    ioi_count += 1
                    if s.chord_spacing_index is not np.nan:
                        self._chord_spacing_index_avg += Decimal(s.chord_spacing_index) * s.duration
                        valid_durations_for_spacing += s.duration
                    self._ps_avg += s.ps
                    if self._ps_max < s.ps:
                        self._ps_max = s.ps
                    if self._ps_min > s.ps:
                        self._ps_min = s.ps
                    if self._pitch_lowest > s.pseg[0].p:
                        self._pitch_lowest = s.pseg[0].p
                    if self._pitch_highest < s.pseg[len(s.pseg) - 1].p:
                        self._pitch_highest = s.pseg[len(s.pseg) - 1].p
                    for v in range(self._num_voices):
                        if len(s.psegs[v]) > 0:
                            # print(v, len(self._pitch_lowest_voices), len(s.psegs))
                            if self._pitch_lowest_voices[v] > s.psegs[v][0].p:
                                self._pitch_lowest_voices[v] = s.psegs[v][0].p
                            if self._pitch_highest_voices[v] < s.psegs[v][len(s.psegs[v]) - 1].p:
                                self._pitch_highest_voices[v] = s.psegs[v][len(s.psegs[v]) - 1].p
                if s.uns is not None:
                    self._ins_avg += s.ins
                    self._lns_avg += s.lns
                    self._median_trajectory_avg += s.median_trajectory
                    self._uns_avg += s.uns
                    if self._ins_max < s.ins:
                        self._ins_max = s.ins
                    if self._ins_min > s.ins:
                        self._ins_min = s.ins
                    if self._lns_max < s.lns:
                        self._lns_max = s.lns
                    if self._lns_min > s.lns:
                        self._lns_min = s.lns
                    if self._median_trajectory_max < s.median_trajectory:
                        self._median_trajectory_max = s.median_trajectory
                    if self._median_trajectory_max > s.median_trajectory:
                        self._median_trajectory_max = s.median_trajectory
                    if self._uns_max < s.uns:
                        self._uns_max = s.uns
                    if self._uns_min > s.uns:
                        self._uns_min = s.uns
                if self._max_pitch_count_with_duplicates < s.pset_cardinality:
                    self._max_pitch_count_with_duplicates = s.pset_cardinality

            # Calculate pitch frequency
            for i in range(0, len(self._slices)):
                for v in range(len(self._slices[i].psets)):
                    for p in self._slices[i].psets[v]:
                        if p.p not in self._pitch_duration_voices[v].keys():
                            self._pitch_duration_voices[v][p.p] = self._slices[i].duration
                        else:
                            self._pitch_duration_voices[v][p.p] += self._slices[i].duration
                        if p.p not in self._pitch_frequency_voices[v].keys():
                            self._pitch_frequency_voices[v][p.p] = 1
                        elif len(self._slices[i-1].psets) == 0 or p not in self._slices[i-1].psets[v]:
                            self._pitch_frequency_voices[v][p.p] += 1
                    for pc in self._slices[i].pcsets[v]:
                        if pc.pc not in self._pc_duration_voices[v].keys():
                            self._pc_duration_voices[v][pc.pc] = self._slices[i].duration
                        else:
                            self._pc_duration_voices[v][pc.pc] += self._slices[i].duration
                        if pc.pc not in self._pc_frequency_voices[v].keys():
                            self._pc_frequency_voices[v][pc.pc] = 1
                        elif len(self._slices[i - 1].pcsets) == 0 or pc not in self._slices[i - 1].pcsets[v]:
                            self._pc_frequency_voices[v][pc.pc] += 1
                for p in self._slices[i].pset:
                    if p.p not in self._pitch_duration.keys():
                        self._pitch_duration[p.p] = self._slices[i].duration
                    else:
                        self._pitch_duration[p.p] += self._slices[i].duration
                    if p.p not in self._pitch_frequency.keys():
                        self._pitch_frequency[p.p] = 1
                    elif p not in self._slices[i-1].pset:
                        self._pitch_frequency[p.p] += 1
                for pc in self._slices[i].pcset:
                    if pc.pc not in self._pc_duration.keys():
                        self._pc_duration[pc.pc] = self._slices[i].duration
                    else:
                        self._pc_duration[pc.pc] += self._slices[i].duration
                    if pc.pc not in self._pc_frequency.keys():
                        self._pc_frequency[pc.pc] = 1
                    elif pc not in self._slices[i-1].pcset:
                        self._pc_frequency[pc.pc] += 1

            # Calculate cseg frequency
            for s in self._slices:
                cseg = s.get_chord_spacing_contour_string()
                pset = s.get_pset_string()
                psc = str(s.ipseg)
                if cseg not in self._chord_spacing_contour_frequency:
                    self._chord_spacing_contour_frequency[cseg] = 1
                    self._chord_spacing_contour_duration[cseg] = s.duration
                else:
                    self._chord_spacing_contour_frequency[cseg] += 1
                    self._chord_spacing_contour_duration[cseg] += s.duration
                if pset not in self._pset_frequency:
                    self._pset_frequency[pset] = 1
                    self._pset_duration[pset] = s.duration
                else:
                    self._pset_frequency[pset] += 1
                    self._pset_duration[pset] += s.duration
                if psc not in self._psc_frequency:
                    self._psc_frequency[psc] = 1
                    self._psc_duration[psc] = s.duration
                else:
                    self._psc_frequency[psc] += 1
                    self._psc_duration[psc] += s.duration
                if s.sc_name not in self._pcsc_frequency:
                    self._pcsc_frequency[s.sc_name] = 1
                    self._pcsc_duration[s.sc_name] = s.duration
                else:
                    self._pcsc_frequency[s.sc_name] += 1
                    self._pcsc_duration[s.sc_name] += s.duration

        # Finalize average calculation
        self._ps_avg /= len(self._slices)
        self._pset_card_avg = float(self._pset_card_avg / self._duration)
        self._pitch_count_with_duplicates_avg = float(self._pitch_count_with_duplicates_avg / self.duration)
        self._ioi_avg_in_seconds /= ioi_count
        self._duration_avg = self.duration / len(self._slices)
    
        if valid_durations_for_spacing != 0:
            self._chord_spacing_index_avg = float(self._chord_spacing_index_avg / valid_durations_for_spacing)
        else:
            self._chord_spacing_index_avg = 0.5

        non_null = self._get_non_null()
        if non_null != 0:
            self._ins_avg /= non_null
            self._lns_avg /= non_null
            self._median_trajectory_avg /= non_null
            self._uns_avg /= non_null

    def _get_non_null(self):
        """
        Gets the number of slices that do not contain None for LNS and UNS
        :return: The number of slices
        """
        counter = 0
        for s in self._slices:
            if s.uns is not None:
                counter += 1
        return counter
