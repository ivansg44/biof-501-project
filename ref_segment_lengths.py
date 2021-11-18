"""Generates json file containing reference segment lengths."""

import json

from Bio import SeqIO

# All fasta files in the reference directory
FASTA_FILES = snakemake.input
# SeqIO objects for each fasta file
SEQ_RECORDS = [next(SeqIO.parse(e, "fasta")) for e in FASTA_FILES]

ref_segment_lengths = {}
for seq_record in SEQ_RECORDS:
    # Parse fasta header to determine which segment we are iterating
    # over.
    header_first_split = seq_record.description.split("|")
    header_second_split = [e.split(":") for e in header_first_split]
    parsed_header_dict = {k: v for [k, v] in header_second_split}

    segment = parsed_header_dict["Segment"]
    ref_segment_lengths[segment] = len(seq_record)

with open(snakemake.output[0], "w") as fp:
    json.dump(ref_segment_lengths, fp)
