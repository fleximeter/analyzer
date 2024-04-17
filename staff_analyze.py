"""
File: staff_analyze.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functions for analyzing individual staves in a score.
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


import music21
import pandas as pd


def parser(score):
    """
    Parses a score and gets the individual note, rest, and chord items. Requires that each staff
    not have separate Voice objects.
    :param score: The score to parse
    :return: A list of staves with note, rest, and chord objects and associated IOI information
    """
    tempo_list = []
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
                    for k, item in enumerate(measure):
                        if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(item) == music21.chord.Chord:
                            stream_collection[staff_idx].append(item)
                        elif type(item) == music21.tempo.MetronomeMark:
                            tempo_list.append((measure.number, item))
            staff_idx += 1
    
    if len(tempo_list) == 0:
        tempo_list.append((0, music21.tempo.MetronomeMark(number=60, referent=1)))

    # Get a better representation of the staves, combining notes that belong together, and assigning duration in seconds
    for i, staff in enumerate(stream_collection):
        staff_collection.append([])
        currentTempoIdx = 0
        for j, item in enumerate(staff):
            # Update the current tempo if necessary
            if len(tempo_list) > currentTempoIdx + 1 and item.measureNumber >= tempo_list[currentTempoIdx + 1][0]:
                currentTempoIdx += 1

            # A simplified representation of the note, with additional information
            note = {
                "measureNumber": item.measureNumber,
                "quarterLength": item.duration.quarterLength,
                "isRest": True if type(item) == music21.note.Rest else False,
                "duration": tempo_list[currentTempoIdx][1].durationToSeconds(item.duration),
                "pitches": item.pitches
            }

            # Continue a tie if necessary; otherwise just add the note
            if item.tie is not None and item.tie.type in ["continue", "stop"]:
                if len(staff_collection[i][-1]["pitches"]) == len(item.pitches):
                    # TEMPORARY HACK
                    valid = True 
                    for pitch in item.pitches:
                        for pitch2 in staff_collection[i][-1]["pitches"]:
                            if pitch.ps == pitch2.ps:
                                valid = True
                    if not valid:
                        raise Exception(f"Bad tie in measure {item.measureNumber}, voice {i}")
                    else:
                        staff_collection[i][-1]["duration"] += note["duration"]
                        staff_collection[i][-1]["quarterLength"] += note["quarterLength"]
            else:
                staff_collection[i].append(note)
        
        # Calculate IOI
        k = None # the note we're calculating IOI for
        ioi = 0
        for j, item in enumerate(staff_collection[i]):
            if item["isRest"] and k:
                ioi += item["duration"]
            elif not item["isRest"]:
                if k is not None:
                    staff_collection[i][k]["ioi"] = ioi
                k = j
                ioi = item["duration"]
    
    return staff_collection


def write_analysis_to_file(path, staff_collection):
    """
    Writes the analysis to a file
    :param path: The file path
    :param staff_collection: The staff collection from the analysis
    """
    dfs = []
    for staff in staff_collection:
        num_pitches = 0
        for item in staff:
            num_pitches = max(num_pitches, len(item["pitches"]))
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
            df["measure_no"].append(item["measureNumber"])
            df["quarter_length"].append(item["quarterLength"])
            df["is_rest"].append(item["isRest"])
            if "ioi" in item:
                df["ioi"].append(item["ioi"])
            else:
                df["ioi"].append(None)
            df["duration"].append(item["duration"])
            for j in range(num_pitches):
                if j < len(item["pitches"]):
                    df[f"pitch{j+1}"].append(item["pitches"][j].ps)        
                else:
                    df[f"pitch{j+1}"].append(None)
        dfs.append(pd.DataFrame(df))
    
    with pd.ExcelWriter(path=path, engine="xlsxwriter") as writer:
        for i, df in enumerate(dfs):
            df.to_excel(writer, sheet_name=f"staff_{i+1}", index=False)
    