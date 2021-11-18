"""Generates visualization file for inputted segment.

This file contains multiple histograms corresponding to different
mutation types across difference strains.
"""

import json

from Bio import AlignIO
import matplotlib.pyplot as plt

with open("reference-segment-lengths.json") as fp:
    REF_SEGMENT_LENGTHS = json.load(fp)

SEGMENT = snakemake.wildcards.segment
SEGMENT_LENGTH = REF_SEGMENT_LENGTHS[SEGMENT]

# matplotlib histogram bins. See https://bit.ly/3DriXDl for details.
BINS = range(1, SEGMENT_LENGTH + 101, 100)

# All inputted pairwise alignment output files
NEEDLE_FILES = [open(e) for e in snakemake.input]
# AlignIO objects corresponding to alignment files
ALIGNMENTS = [AlignIO.read(e, "emboss") for e in NEEDLE_FILES]

fig, axs = plt.subplots(
    # One row for each strain; which corresponds to the number of
    # alignments.
    len(ALIGNMENTS),
    # One col for each mutation types
    3,
    # Minimize padding
    tight_layout=True,
    # Consistent y axis scale in each col
    sharey="col",
    # Fig size correlates with the number of bins, which correlates
    # with inputted segment length.
    figsize=(len(BINS) * 0.5, 6)
)

# Title at top of visualization file
fig.suptitle("Segment %s" % SEGMENT)

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

    # Lists of reference genome positions where mutations occur in the
    # alternate genome.
    snps = []
    insertions = []
    deletions = []
    # Insertions in alternate genome lead to frame shift
    frame_shift = 0

    for j in range(start_index, end_index+1):
        ref_index = j - frame_shift
        if ref_seq[j] == "-":
            insertions.append(ref_index)
            frame_shift += 1
        elif alt_seq[j] == "-":
            deletions.append(ref_index)
        elif ref_seq[j] != alt_seq[j]:
            snps.append(ref_index)

    # Histogram for each type of mutation, in row corresponding to this
    # alignment.
    axs[i, 0].hist(snps, bins=BINS, rwidth=0.7, color="black")
    axs[i, 1].hist(insertions, bins=BINS, rwidth=0.7, color="#4daf4a")
    axs[i, 2].hist(deletions, bins=BINS, rwidth=0.7, color="#e41a1c")

    # Put a label on the left of this row indicating which alternate
    # strain was used.
    alt_id = alignment[1].id
    plt.setp(axs[i], ylabel=alt_id)

# Set titles at the top of each column to indicate what mutations are
# visualized.
axs[0, 0].set_title("SNPs")
axs[0, 1].set_title("Insertions")
axs[0, 2].set_title("Deletions")

# No histogram should have negative values. If we do not set this, and
# there are no mutations in a column, the histogram y-axis will center
# at 0 with negative values below it.
plt.setp(axs, ylim=0)

plt.savefig(snakemake.output[0])
