"""Generates visualization file for inputted segment.

This file contains multiple heatmaps corresponding to different
mutation types across different strains.
"""

import json
from math import floor
from os import path

from Bio import AlignIO
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np

with open(path.join("assets", "reference-segment-lengths.json")) as fp:
    REF_SEGMENT_LENGTHS = json.load(fp)

SEGMENT = snakemake.wildcards.segment
SEGMENT_LENGTH = REF_SEGMENT_LENGTHS[SEGMENT]

# We will bin mutations by nt pos along x axis
BIN_EDGES = range(1, SEGMENT_LENGTH + 101, 100)
BIN_LABELS = ["%s - %s" % (e, e+100) for e in BIN_EDGES[:-1]]

# All inputted pairwise alignment output files
NEEDLE_FILES = [open(e) for e in snakemake.input]
# AlignIO objects corresponding to alignment files
ALIGNMENTS = [AlignIO.read(e, "emboss") for e in NEEDLE_FILES]
# Alternate strains from each alignment
STRAINS = [e[1].id for e in ALIGNMENTS]

with open(path.join("assets", "pandemic-strains.json")) as fp:
    PANDEMIC_STRAINS = json.load(fp)
# Indexes in STRAINS corresponding to pandemic strains
PANDEMIC_INDICES = [i for i, e in enumerate(STRAINS) if e in PANDEMIC_STRAINS]

fig, axs = plt.subplots(
    # One row for each mutation type
    3,
    # Only one heatmap per row
    1,
    # Fig size correlates with the number of bins and strains. The
    # latter correlates with number of alignments.
    figsize=(len(BIN_EDGES) * 0.5, len(ALIGNMENTS))
)

# Title at top of visualization file
fig.suptitle("Segment %s" % SEGMENT)

# https://stackoverflow.com/a/45161551/11472358
fig.tight_layout(h_pad=8, rect=[0.06, 0.06, 1, 0.95])

# Iterate over subplots to implement stylistic changes prior to
# generating heatmap data.
for i in range(3):
    axs[i].set_xticks(np.arange(len(BIN_LABELS)))
    axs[i].set_yticks(np.arange(len(ALIGNMENTS)))
    axs[i].set_xticklabels(BIN_LABELS, rotation=270)
    axs[i].set_yticklabels(STRAINS)

    # Color pandemic tick labels red
    for j in PANDEMIC_INDICES:
        axs[i].get_yticklabels()[j].set_color("red")

    if i == 0:
        axs[i].set_title("SNPs")
    elif i == 1:
        axs[i].set_title("Insertions")
    elif i == 2:
        axs[i].set_title("Deletions")

# Matrices that will count the number of mutations in each bin relative
# to reference genome position.
snp_matrix = \
    np.zeros(shape=(len(ALIGNMENTS), len(BIN_LABELS)), dtype=int)
insertion_matrix = \
    np.zeros(shape=(len(ALIGNMENTS), len(BIN_LABELS)), dtype=int)
deletion_matrix = \
    np.zeros(shape=(len(ALIGNMENTS), len(BIN_LABELS)), dtype=int)

for i, alignment in enumerate(ALIGNMENTS):
    # The alignment is ultimately derived from a reference and
    # alternate fasta file.
    ref_seq, alt_seq = alignment[0].seq, alignment[1].seq

    # All mutation positions are marked relative to the reference
    # genome, so we must determine where the reference genome begins
    # and ends. The alternate genome may have extra bases at the
    # beginning or end.
    start_index = 0
    for j, base in enumerate(ref_seq):
        if base != "-":
            start_index = j
            break
    end_index = -1
    for j, base in reversed(list(enumerate(ref_seq))):
        if base != "-":
            end_index = j
            break

    # Insertions in alternate genome lead to frame shift
    frame_shift = 0

    for j in range(start_index, end_index+1):
        # Index relative to reference genome position
        ref_index = j - start_index - frame_shift
        ref_index_bin = int(floor(ref_index/100))

        if ref_seq[j] == "-":
            insertion_matrix[i, ref_index_bin] += 1
            frame_shift += 1
        elif alt_seq[j] == "-":
            deletion_matrix[i, ref_index_bin] += 1
        elif ref_seq[j] != alt_seq[j]:
            snp_matrix[i, ref_index_bin] += 1
        j += 1

# Label heatmap cells with the actual number of mutations
matrices = [snp_matrix, insertion_matrix, deletion_matrix]
for i, matrix in enumerate(matrices):
    for j, x in enumerate(matrix):
        for k, y in enumerate(x):
            axs[i].text(k, j, matrix[j, k],
                        ha="center", va="center", color="w",
                        path_effects=[
                            path_effects.Stroke(linewidth=1,
                                                foreground='black'),
                            path_effects.Normal()
                        ])

# Specify custom color schemes
axs[0].imshow(snp_matrix, cmap=plt.get_cmap("Blues"))
axs[1].imshow(insertion_matrix, cmap=plt.get_cmap("Greens"))
axs[2].imshow(deletion_matrix, cmap=plt.get_cmap("Reds"))

plt.savefig(snakemake.output[0])
