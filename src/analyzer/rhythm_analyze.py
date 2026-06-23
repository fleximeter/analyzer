"""
File: rhythm_analyze.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functions for analyzing rhythm in individual staves in a score.
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


from fractions import Fraction
import music21
import pandas as pd


class CollapsedNote:
    """
    Represents a note, chord, or rest that has been collapsed from tied items, including across barlines
    """
    def __init__(self, *args, **kwargs):
        self._is_rest = kwargs["is_rest"] if "is_rest" in kwargs else False
        self._measure_number = kwargs["measure_number"] if "measure_number" in kwargs else 0
        self._quarter_length = kwargs["quarter_length"] if "quarter_length" in kwargs else 0
        self._duration = kwargs["duration"] if "duration" in kwargs else 0
        self._pitches = kwargs["pitches"] if "pitches" in kwargs else None
        self._start_time = kwargs["start_time"] if "start_time" in kwargs else 0
        self._ioi = kwargs["ioi"] if "ioi" in kwargs else None
    
    @property
    def is_rest(self):
        """
        Whether or not the current CollapsedNote represents a rest
        """
        return self._is_rest

    @property
    def measure_number(self):
        """
        The measure number of the CollapsedNote
        """
        return self._measure_number

    @property
    def quarter_length(self):
        """
        The length of the CollapsedNote, in quarter notes
        """
        return self._quarter_length

    @quarter_length.setter
    def quarter_length(self, value):
        """
        The length of the CollapsedNote, in quarter notes
        """
        self._quarter_length = value

    @property
    def duration(self):
        """
        The duration of the CollapsedNote
        """
        return self._duration

    @duration.setter
    def duration(self, value):
        """
        The duration of the CollapsedNote
        """
        self._duration = value

    @property
    def pitches(self):
        """
        The pitches of the CollapsedNote
        """
        return self._pitches

    @property
    def start_time(self):
        """
        The start time of the CollapsedNote
        """
        return self._start_time

    @property
    def ioi(self):
        """
        The IOI of the CollapsedNote, in seconds
        """
        return self._ioi

    @ioi.setter
    def ioi(self, value):
        """
        The IOI of the CollapsedNote, in seconds
        """
        self._ioi = value


def find_duration_successsion(ioi_analysis, succession, include_rests=False):
    """
    Finds a rhythm duration succession in an IOI analysis (such as (1/2, 1/4, 1/3, 2)).
    These are hard quarter lengths. They require specific note values and do not flag 
    augmentation and diminution.
    
    The algorithm records the measure numbers in which these successions start for each voice.
    If a succession starts multiple times in the same measure, the measure will be recorded
    multiple times. The results are returned in a dictionary, with separate counts for each staff.

    :param ioi_analysis: The IOI analysis
    :param succession: A tuple or list of numbers. It is best to use integers or the Fraction type.
    :param include_rests: Whether or not to count rests as part of the duration succession
    :return: A dictionary of rhythm successions
    """
    voice_dict = {}
    for v, voice in enumerate(ioi_analysis):
        counter = 0
        measures = []
        for i in range(len(voice) - len(succession) + 1):
            succession_found = True
            for j in range(0, len(succession)):
                # If the current note is a rest, the succession is invalid
                if voice[i+j].is_rest and not include_rests:
                    succession_found = False
                    break
                # If the current ratio is wrong, the succession is invalid
                elif voice[i+j].quarter_length != succession[j]:
                    succession_found = False
                    break
            if succession_found:
                counter += 1
                measures.append(voice[i].measure_number)
        voice_dict[f"voice_{v+1}"] = {"num": counter, "measures_start": measures}
    return voice_dict


def find_rhythm_ratio_successsion(ioi_analysis, succession, include_rests=False):
    """
    Finds a rhythm ratio succession in an IOI analysis (such as (2, 1, 3, 1, 2, 1)).
    The ratios are simply ratios of duration. They do not require specific note values,
    For example, the succession given could be realized as half, quarter, dotted half, etc.
    or as quarter, eighth, dotted quarter, etc. Both would be flagged using this algorithm.

    The algorithm records the measure numbers in which these successions start for each voice.
    If a succession starts multiple times in the same measure, the measure will be recorded
    multiple times. The results are returned in a dictionary, with separate counts for each staff.

    :param ioi_analysis: The IOI analysis
    :param succession: A tuple or list of numbers. It is best to use integers or the Fraction type.
    :param include_rests: Whether or not to count rests as part of the rhythm ratio succession
    :return: A dictionary of rhythm ratio successions
    """
    voice_dict = {}
    for v, voice in enumerate(ioi_analysis):
        counter = 0
        measures = []
        for i in range(len(voice) - len(succession) + 1):
            succession_found = True

            # If the current note is a rest, the succession is invalid
            if voice[i].is_rest and not include_rests:
                succession_found = False
            else:
                for j in range(1, len(succession)):
                    # If the current note is a rest, the succession is invalid
                    if voice[i+j].is_rest and not include_rests:
                        succession_found = False
                        break
                    # If the current ratio is wrong, the succession is invalid
                    elif voice[i+j].quarter_length != voice[i+j-1].quarter_length * Fraction(succession[j], succession[j-1]):
                        succession_found = False
                        break
            if succession_found:
                counter += 1
                measures.append(voice[i].measure_number)
        voice_dict[f"voice_{v+1}"] = {"num": counter, "measures_start": measures}
    return voice_dict


def part_ioi_analyzer(score, tempo_list=None):
    """
    Parses a score and gets the individual note, rest, and chord items. Requires that each staff
    not have separate Voice objects.
    :param score: The score to parse
    :return: A list of staves with note, rest, and chord objects and associated IOI information
    """
    if tempo_list is None:
        tempo_list_from_score = []
    stream_collection = []
    staff_collection = []

    score = music21.converter.parse(score)

    # Get the notes in a nice easy structure
    staff_idx = 0
    for i, staff in enumerate(score):
        if type(staff) in [music21.stream.Part, music21.stream.PartStaff]:
            stream_collection.append([])
            for j, measure in enumerate(staff):                
                if type(measure) == music21.stream.Measure:
                    for current_note, item in enumerate(measure):
                        if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(item) == music21.chord.Chord:
                            stream_collection[staff_idx].append(item)
                        elif type(item) == music21.tempo.MetronomeMark:
                            tempo_list_from_score.append((measure.number, item))
            staff_idx += 1
    
    # What if no tempo is provided?
    if tempo_list is None and len(tempo_list_from_score) == 0:
        tempo_list_from_score.append((0, music21.tempo.MetronomeMark(number=60, referent=1)))
    if tempo_list is not None:
        tempo_list_from_score = [(tempo, music21.tempo.MetronomeMark(number=tempo_list[tempo], referent=1)) for tempo in tempo_list]

    # Get a better representation of the staves, combining notes that belong together, and assigning duration in seconds
    for i, staff in enumerate(stream_collection):
        staff_collection.append([])
        currentTempoIdx = 0
        for j, item in enumerate(staff):
            start_time = 0

            # Update the current tempo if necessary
            if len(tempo_list_from_score) > currentTempoIdx + 1 and item.measureNumber >= tempo_list_from_score[currentTempoIdx + 1][0]:
                currentTempoIdx += 1

            # A simplified representation of the note, with additional information
            note = CollapsedNote(
                measure_number=item.measureNumber,
                quarter_length=item.duration.quarterLength,
                is_rest=True if type(item) == music21.note.Rest else False,
                duration=tempo_list_from_score[currentTempoIdx][1].durationToSeconds(item.duration),
                pitches=item.pitches,
                start_time=start_time
            )

            start_time += note.duration

            # Continue a tie if necessary; otherwise just add the note
            if item.tie is not None and item.tie.type in ["continue", "stop"]:
                if len(staff_collection[i][-1].pitches) == len(item.pitches):
                    valid = True 
                    for pitch in item.pitches:
                        for pitch2 in staff_collection[i][-1].pitches:
                            if pitch.ps == pitch2.ps:
                                valid = True
                    if not valid:
                        raise Exception(f"Bad tie in measure {item.measureNumber}, voice {i}")
                    else:
                        staff_collection[i][-1].duration += note.duration
                        staff_collection[i][-1].quarter_length += note.quarter_length
            else:
                staff_collection[i].append(note)
        
        # Calculate IOI
        current_note = None # the note we're calculating IOI for
        current_ioi = 0
        for j, item in enumerate(staff_collection[i]):
            if item.is_rest:
                current_ioi += item.duration
            else:
                # update ioi
                if current_note is not None:
                    current_note.ioi = current_ioi
                # reset current note
                current_note = item
                current_ioi = item.duration
        # catch the last note
        if current_note is not None:
            current_note.ioi = current_ioi        
    
    return staff_collection


def write_ioi_analysis_to_file(path, staff_collection):
    """
    Writes the analysis to a file
    :param path: The file path
    :param staff_collection: The staff collection from the analysis
    """
    dfs = []
    for staff in staff_collection:
        num_pitches = 0
        for item in staff:
            num_pitches = max(num_pitches, len(item.pitches))
        df = {
            "measure_no": [],
            "quarter_length": [],
            "is_rest": [],
            "ioi": [],
            "duration": []
        }
        for j in range(num_pitches):
            df[f"pitch{j+1}"] = []
        for item in staff:
            df["measure_no"].append(item.measure_number)
            df["quarter_length"].append(item.quarter_length)
            df["is_rest"].append(item.is_rest)
            df["ioi"].append(item.ioi)
            df["duration"].append(item.duration)
            for j in range(num_pitches):
                if j < len(item.pitches):
                    df[f"pitch{j+1}"].append(item.pitches[j].ps)        
                else:
                    df[f"pitch{j+1}"].append(None)
        dfs.append(pd.DataFrame(df))
    
    with pd.ExcelWriter(path=path, engine="xlsxwriter") as writer:
        for i, df in enumerate(dfs):
            df.to_excel(writer, sheet_name=f"staff_{i+1}", index=False)
    