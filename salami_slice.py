"""
File: salami_slice.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains the SalamiSlice class for salami slicing with music21.
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

from pctheory import cseg, pitch
from decimal import Decimal
import numpy as np


def sort_pitch_name_list(pitch_name_list):
    """
    Sorts a list of pitch names
    :param pitch_name_list: A list of pitch names
    :return: None
    """
    lmap = {'C': 1, 'D': 4, 'E': 7, 'F': 10, 'G': 13, 'A': 16, 'B': 19}
    pitch_name_list2 = []
    for p in pitch_name_list:
        p1 = int(p[-1]) * 21 + lmap[p[0]]
        if len(p) > 2:
            if p[1] == '#':
                p1 += 1
            elif p[1] == '-':
                p1 -= 1
        pitch_name_list2.append((p, p1))
    return [p[0] for p in sorted(pitch_name_list2, key=lambda x: x[1])]


class SalamiSlice:
    def __init__(self, tempo=1, quarter_duration=1, measure=None, num_voices=1):
        """
        Creates a SalamiSlice
        :param duration: The duration of the slice, in seconds
        :param quarter_duration: The duration of the slice, in quarter notes
        :param measure: The measure number
        :param aslice: An existing slice if using copy constructor functionality
        """
        self._core = False                      # Whether or not the chord is a core harmony
        self._chord_spacing_contour = None      # The contour of the pset
        self._chord_spacing_index = 0           # The spread measure of the chord
        self._derived_core = False              # Whether or not the chord is a derived core harmony
        self._derived_core_associations = None  # Derived core associations, if any
        self._duration = 0                      # The duration of the slice in seconds
        self._ipseg = []                        # The ipseg of the slice
        self._ioi_in_seconds = None                        # The IOI
        self._measure = measure                 # The measure number in which the slice begins
        self._num_voices = num_voices           # The number of voices
        self._pset_cardinality = 0              # The number of distinct pitches present (excluding duplicates)
        self._pitch_count_with_duplicates = 0   # The number of pitches present (including duplicates)
        self._pcset_cardinality = 0             # The number of distinct pitch-classes present (excluding duplicates)
        self._pcseg = None                      # The pcseg
        self._pcsegs = [[] for i in range(num_voices)]      # The pcsegs by voice
        self._pcset = None                      # The pcset
        self._pcsets = [set() for i in range(num_voices)]   # The pcsets by voice
        self._pitchseg = []                     # The pitch seg (contains duplicates)
        self._pitchsegs = [[] for i in range(num_voices)]   # The pitch segs by voice
        self._pitch_name_list = []                     # A list of pitch names
        self._pitch_name_lists = [[] for i in range(num_voices)]   # The pitch_name_lists by voice
        self._pset = None                       # The pset
        self._psets = [set() for i in range(num_voices)]    # The psets by voice
        self._pseg = None                       # The pseg
        self._psegs = []                        # The psegs by voice
        self._quarter_duration = quarter_duration           # The duration in quarters
        self._sc_name = None                    # The set-class name of the pcset
        self._sc_name_carter = None             # The Carter set-class name of the pcset
        self._tempo = tempo                     # The tempo

        # Calculations after Burt (2012)
        self._ins = None                        # The INS of the slice
        self._lns = None                        # The LNS of the slice
        self._lower_bound = None                # The lower bound of the slice.
        self._median_trajectory = None          # The MT of the slice
        
        # The NS of the slice. If the lower and upper bounds are defined, but LNS and UNS are not,
        # the NS represents the entire pitch area encompassed by the piece. Otherwise it is None.
        self._ns = None
        
        self._ps = None                         # The PS of the slice
        self._start_position = None             # The start position in quarter notes relative to the current measure
        self._time_signature = None             # The time signature of the slice
        self._uns = None                        # The UNS of the slice
        self._upper_bound = None                # The upper bound of the slice.

    @property
    def chord_spacing_contour(self):
        """
        The contour of the pseg
        :return: The contour
        """
        return self._chord_spacing_contour

    @property
    def chord_spacing_index(self):
        """
        The pset spacing index of the SalamiSlice
        :return: The pset spacing index
        """
        return self._chord_spacing_index

    @property
    def core(self):
        """
        Whether or not the chord is a core harmony
        :return: True or False
        """
        return self._core

    @property
    def derived_core(self):
        """
        Whether or not the chord is a derived core harmony
        :return: True or False
        """
        return self._derived_core

    @property
    def derived_core_associations(self):
        """
        Derived core associations, if any
        :return: Derived core associations, if any
        """
        return self._derived_core_associations

    @property
    def duration(self):
        """
        The duration
        :return: The duration
        """
        return self._duration

    @duration.setter
    def duration(self, value):
        """
        The duration
        :param value: The new duration
        :return: None
        """
        self._duration = value

    @property
    def ins(self):
        """
        The internal negative space (INS)
        :return: The internal negative space (INS)
        """
        return self._ins
    
    @property
    def ioi_in_seconds(self):
        """
        The interonset interval (IOI) in *seconds*
        :return: The interonset interval (IOI) in *seconds*
        """
        return self._ioi_in_seconds

    @property
    def ipseg(self):
        """
        The ipseg (ordered interval succession between adjacent pitches from low to high)
        :return: The ipseg
        """
        return self._ipseg

    @property
    def lns(self):
        """
        The lower negative space (LNS)
        :return: The lower negative space (LNS)
        """
        return self._lns

    @property
    def lower_bound(self):
        """
        The lower bound
        :return: The lower bound
        """
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        """
        The lower bound
        :param value: The new lower bound
        :return: None
        """
        self._lower_bound = value

    @property
    def measure(self):
        """
        The measure number
        :return: The measure number
        """
        return self._measure

    @property
    def median_trajectory(self):
        """
        The median trajectory (MT)
        :return: The median trajectory (MT)
        """
        return self._median_trajectory

    @property
    def ns(self):
        """
        The negative space (NS)
        :return: The negative space (NS)
        """
        return self._ns

    @property
    def pset_cardinality(self):
        """
        The pitch cardinality of the SalamiSlice (excludes duplicates)
        :return: The pitch cardinality
        """
        return self._pset_cardinality

    @property
    def pitch_count_with_duplicates(self):
        """
        The pitch count of the SalamiSlice (contains duplicates)
        :return: The pitch count
        """
        return self._pitch_count_with_duplicates

    @property
    def pcset_cardinality(self):
        """
        The pitch-class set cardinality of the SalamiSlice (excludes duplicates)
        :return:
        """
        return self._pcset_cardinality

    @property
    def ps(self):
        """
        The positive space (PS)
        :return: The positive space (PS)
        """
        return self._ps

    @property
    def pcseg(self):
        """
        The pcseg of the SalamiSlice
        :return: The pcseg
        """
        return self._pcseg

    @property
    def pcset(self):
        """
        The pcset of the SalamiSlice
        :return: The pcset
        """
        return self._pcset

    @property
    def pcsegs(self):
        """
        The pcseg of the SalamiSlice
        :return: The pcseg
        """
        return self._pcsegs

    @property
    def pcsets(self):
        """
        The pcset of the SalamiSlice
        :return: The pcset
        """
        return self._pcsets

    @property
    def pitchseg(self):
        """
        The pitchseg of the SalamiSlice
        :return: The pitchseg
        """
        return self._pitchseg

    @property
    def pitch_name_list(self):
        """
        A list of pitch names
        :return: A list of pitch names
        """
        return self._pitch_name_list

    @property
    def pitchsegs(self):
        """
        The pitchseg of the SalamiSlice
        :return: The pitchseg
        """
        return self._pitchsegs

    @property
    def pitch_name_lists(self):
        """
        A list of pitch names
        :return: A list of pitch names
        """
        return self._pitch_name_lists

    @property
    def pseg(self):
        """
        The pseg of the SalamiSlice (contains duplicates)
        :return: The pseg
        """
        return self._pseg

    @property
    def pset(self):
        """
        The pset of the SalamiSlice
        :return: The pset
        """
        return self._pset

    @property
    def psegs(self):
        """
        The pseg of the SalamiSlice (contains duplicates)
        :return: The pseg
        """
        return self._psegs

    @property
    def psets(self):
        """
        The pset of the SalamiSlice
        :return: The pset
        """
        return self._psets

    @property
    def quarter_duration(self):
        """
        The duration in quarter notes
        :return: The duration in quarter notes
        """
        return self._quarter_duration

    @quarter_duration.setter
    def quarter_duration(self, value):
        """
        The duration in quarter notes
        :param value: The new duration
        :return: None
        """
        self._quarter_duration = value

    @property
    def sc_name(self):
        """
        The set-class name of the SalamiSlice
        :return: The set-class name
        """
        return self._sc_name

    @property
    def sc_name_carter(self):
        """
        The Carter name of the SalamiSlice
        :return: The Carter name
        """
        return self._sc_name_carter

    @property
    def start_position(self):
        """
        The start position in the measure
        :return: The start position in the measure
        """
        return self._start_position

    @start_position.setter
    def start_position(self, value):
        """
        The start position in the measure
        :param value: The new start position
        :return: None
        """
        self._start_position = value

    @property
    def time_signature(self):
        """
        The time signature (music21.meter.TimeSignature)
        :return: The time signature
        """
        return self._time_signature

    @time_signature.setter
    def time_signature(self, value):
        """
        The time signature (music21.meter.TimeSignature)
        :param value: The new time signature
        :return: None
        """
        self._time_signature = value

    @property
    def uns(self):
        """
        The upper negative space (UNS)
        :return: The upper negative space (UNS)
        """
        return self._uns

    @property
    def upper_bound(self):
        """
        The upper bound
        :return: The upper bound
        """
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        """
        The upper bound
        :param value: The new upper bound
        :return: None
        """
        self._upper_bound = value

    def add_pitches(self, pitches, pitch_names=None, voice=0):
        """
        Adds pitches to the SalamiSlice
        :param pitches: A collection of pitches to add
        :param pitch_names: A collection of pitch names (as strings) corresponding to the pitch collection
        :param voice: The voice
        """
        # Add each pitch to the chord
        for i in range(len(pitches)):
            self._pitchseg.append(pitches[i])
            self._pitchsegs[voice].append(pitches[i])
            self._pitch_name_list.append(pitch_names[i])
            self._pitch_name_lists[voice].append(pitch_names[i])
            self._pitch_count_with_duplicates += 1

    def copy(self):
        """
        Copies the SalamiSlice
        :return: A copy of the SalamiSlice
        """
        v = SalamiSlice(self._tempo, self._quarter_duration, self._measure, self._num_voices)

    def get_chord_spacing_contour_string(self):
        """
        Gets the chord spacing contour as a string
        :return: The chord spacing contour as a string
        """
        cseg = "<"
        for cp in self._chord_spacing_contour:
            cseg += str(cp) + ", "
        if cseg[-1] == " ":
            cseg = cseg[:-2]
        cseg += ">"
        return cseg

    def get_ipseg_string(self):
        """
        Gets the ipseg as a string
        :return: The ipseg as a string
        """
        ipseg = "<"
        for ip in self._ipseg:
            ipseg += str(ip) + ", "
        if ipseg[-1] == " ":
            ipseg = ipseg[:-2]
        ipseg += ">"
        return ipseg

    def get_pcset_string(self):
        """
        The pcset
        :return: The pcset
        """
        pcset_str = "{"
        sorted_pcset = list(self._pcset)
        sorted_pcset.sort()
        for pc in sorted_pcset:
            pcset_str += str(pc)
        pcset_str += "}"
        return pcset_str

    def get_pset_string(self):
        """
        The pset
        :return: The pset
        """
        pset_str = "{"
        for p in self._pseg:
            pset_str += f"{str(p.p)}, "
        pset_str = pset_str.strip(' ,')
        pset_str += "}"
        return pset_str

    def prepare_for_clean(self):
        """
        Sorts the pitchseg in preparation for slice cleaning
        :return:
        """
        self._pitchseg.sort()

    def run_calculations(self, sc):
        """
        Calculates information about the SalamiSlice. You should combine any SalamiSlices that you want to combine before
        running this method, to avoid making unnecessary computations.
        :param sc: A SetClass object
        :return: None
        """
        self._pset = set()
        self._pcseg = []
        self._pcset = set()
        for p in self._pitchseg:
            self._pset.add(pitch.Pitch(p))
            self._pcset.add(pitch.PitchClass(p))
        for v in range(self._num_voices):
            for p in self._pitchsegs[v]:
                self._psets[v].add(pitch.Pitch(p))
                self._pcsets[v].add(pitch.PitchClass(p))
            self._psegs.append(list(self._psets[v]))
            self._psegs[v].sort()
            for p in self._psegs[v]:
                self._pcsegs[v].append(pitch.PitchClass(p.pc))
        self._pseg = list(self._pset)
        self._pseg.sort()
        self._pitch_name_list = sort_pitch_name_list(self._pitch_name_list)
        for v in range(self._num_voices):
            self._pitchsegs[v].sort()
            self._pitch_name_lists[v] = sort_pitch_name_list(self._pitch_name_lists[v])
        for p in self._pseg:
            self._pcseg.append(pitch.PitchClass(p.pc))
        if len(self._ipseg) > 0:
            self._ipseg.clear()
        for i in range(1, len(self._pseg)):
            self._ipseg.append(self._pseg[i].p - self._pseg[i - 1].p)
        self._chord_spacing_contour = cseg.simplify([ip for ip in self._ipseg])

        # Calculate values
        self._pset_cardinality = len(self._pset)
        self._pcset_cardinality = len(self._pcset)
        self._duration = (Decimal(60) / Decimal(self._tempo)) * (Decimal(self._quarter_duration.numerator) /
                                                                 Decimal(self._quarter_duration.denominator))

        # Calculate pset spacing index
        self._chord_spacing_index = calculate_chord_spacing_index(self._pset)

        # Calculate set theory info
        sc.pcset = self._pcset
        self._sc_name = sc.name_morris
        if 1 < self.pcset_cardinality < 11:
            self._sc_name_carter = sc.name_carter
        else:
            self._sc_name_carter = ""
        if (self._sc_name_carter == "18" or self._sc_name_carter == "23") and self._pcset_cardinality == 4:
            self._core = True
        if self._sc_name_carter == "35" and self._pcset_cardinality == 6:
            self._core = True
        self._derived_core_associations = sc.derived_core
        if self._derived_core_associations is not None:
            self._derived_core = True

    def run_calculations_burt(self):
        """
        Calculates information about the SalamiSlice. You must set the lower and upper bounds before running this
        method. You should also combine any SalamiSlices that you want to combine before running this method,
        to avoid making unnecessary computations.
        :return: None
        """
        # Calculate ps and ins
        self._ps = len(self._pset)
        if self._ps > 0:
            self._ins = self._pseg[-1].p - self._pseg[0].p + 1 - self._ps
        else:
            self._ins = 0

        # Calculate uns, lns, ns, and mt
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._ps is None or self._ps == 0:
                self._lns = None
                self._uns = None
                self._median_trajectory = None
                self._ns = self._upper_bound - self._lower_bound + 1
            else:
                self._lns = self._pseg[0].p - self._lower_bound
                self._uns = self._upper_bound - self._pseg[-1].p
                self._median_trajectory = (self._lns - self._uns) / 2


def calculate_chord_spacing_index(pset):
    """
    Calculates the chord spacing index (CSI) of a pset
    :param pset: The pset
    :return: The CSI
    """
    if len(pset) < 3:
        return np.nan
    else:
        avg = 0
        pseg = sorted(list(pset))
        for p in pseg:
            avg += p.p
        avg /= len(pset)
        return (avg - pseg[0].p) / (pseg[-1].p - pseg[0].p)
