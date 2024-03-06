# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome
from genome_kit import GenomeTrack
from genome_kit import GenomeTrackBuilder
import numpy as np

# Intervals on which our data will be defined
genome = Genome('hg19')
intervals = [
    genome.interval('chr1', '+', 0, 5),
    genome.interval('chr1', '+', 7, 10),
]

# Data to fill those intervals
data = [
    np.array(
        [
            [0, 1],  # For intervals[0]
            [1, 2],
            [2, 4],
            [3, 8],
            [4, 16]
        ],
        np.float16),
    np.array(
        [
            [7, 128],  # For intervals[1]
            [8, 256],
            [9, 512]
        ],
        np.float16),
]

# Build a single-stranded f16 track with 2 values per position.
build = GenomeTrackBuilder('build_track.gtrack', 'f16', 'single_stranded', genome, dim=2)
build.set_default_value(-1)
for interval, data in zip(intervals, data):
    build.set_data(interval, data)
build.finalize()

# Load the track and print the entire range of values.
with GenomeTrack('build_track.gtrack') as track:
    print("As positive strand:\n", track(genome.interval('chr1', '+', 0, 10)))
    print("As negative strand:\n", track(genome.interval('chr1', '-', 0, 10)))
