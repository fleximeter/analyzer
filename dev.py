import pctheory
import music21
from analyzer import salami_slice_analyze
from pathlib import Path

file = Path(__file__).parent / "tests/data/test_carter1.musicxml"
results = salami_slice_analyze.analyze(file)
slices = results[0].slices
print(slices[2].pseg)