# Calibrating against recorded notes

These scripts compare the model with single notes of a real grand piano, note by note
and velocity by velocity, so loss, coupling and felt settings can be checked against
measurements instead of only by ear.

The reference is the Salamander Grand Piano V3 (a Yamaha C5, recorded at 48 kHz with
two AKG C414 about 12 cm above the strings, 16 velocity layers every minor third,
CC-BY 3.0). Download the SFZ+FLAC version (about 740 MB) from
https://freepats.zenvoid.org/Piano/acoustic-grand-piano.html and unpack it anywhere.

```sh
pip install numpy scipy soundfile matplotlib
ninja -C build openpiano-notes
build/tools/openpiano-notes --out renders/notes --seconds 8 \
    --notes 21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108 \
    --vel 40,60,92,124
python3 tools/calibrate/analyze.py --salamander PATH/TO/SalamanderGrandPiano-SFZ+FLAC-V3+20200602 \
    --model renders/notes --out report/calibrate
python3 tools/calibrate/plot.py report/calibrate/metrics.json --vel 92 --out report/calibrate/decay_v92.png
```

`analyze.py` writes `report.md` (one row per note and velocity, real / model) and
`metrics.json`, and prints the mean differences. MIDI velocities 40, 60, 92 and 124 are
the middles of Salamander layers 4, 8, 12 and 16.

What is measured, from the onset:

- **early / late decay** [dB/s]: slope of the harmonic level (the first 20 partials,
  each in a narrow band, summed) over 0.1-1 s and over 2-6 s (earlier for the short
  treble recordings). Narrow bands keep the recordings' hiss out of the late decay.
- **f1 early / late**: the same for the first partial alone.
- **partial decays**: each partial's slope over the first second (in `metrics.json`),
  which shows how losses rise with frequency.
- **brightness**: spectral centroid (in multiples of f1), the share of energy above
  2 kHz, and the level of partials 2-8 against the first, all over 20-200 ms.
- **level**: loudness over the first 200 ms.

Things to keep in mind when reading the numbers:

- The Salamander piano is a 2 m Yamaha C5; the model is built on Steinway D string
  lengths. Matching the decay and the way the tone changes with velocity is meaningful;
  matching every note's exact spectrum is not.
- The microphones sit close above the strings, so the recordings are brighter and the
  bass fundamentals weaker than at the player's seat, where the model listens. The
  velocity-to-velocity changes (how much brighter and louder a note gets when played
  harder) do not depend on the microphones, so they are the better brightness targets.
- One recorded note can have its own quirks (a slightly mistuned unison beats, a
  hot board resonance), so judge trends across neighbouring notes.
