import pctheory
import music21
from analyzer import salami_slice_analyze
from pathlib import Path

file = Path(__file__).parent / "tests/data/test3.musicxml"
results = salami_slice_analyze.analyze(file)
print(len(results[0].slices))

s = music21.converter.parse(file)
for item in s[4][2][1]:
    print(item)