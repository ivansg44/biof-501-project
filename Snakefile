from os import path

SEGMENTS = range(1, 9)
SAMPLES = [
    "A-PuertoRico-8-1934",
    "A-FortMonmouth-1-1947",
    "A-Albany-4836-1950",
    "A-NewJersey-1976",
    "A-USSR-92-1977",
    "A-India-6263-1980",
    "A-Taiwan-1-1986",
    "A-Texas-36-91",
    "A-Beijing-262-1995",
    "A-NewYork-222-2003",
    "A-California-07-2009",
    "A-Alabama-13-2015",
    "A-Alabama-01-2020"
]

rule all:
    input:
        "reference-segment-lengths.json",
        expand("visualizations/{segment}.pdf", segment=SEGMENTS)

rule pairwise_alignment:
    input:
        path.join(
            "reference",
            "segment-{segment}-reference",
            "segment-{segment}-reference.fasta"
        ),
        path.join(
            "samples",
            "segment-{segment}-samples",
            "{sample}-segment-{segment}.fasta"
        )
    output:
        path.join(
            "pairwise-alignments",
            "segment-{segment}-pairwise-alignments",
            "{sample}-segment-{segment}-alignment.needle"
        )
    shell:
        "needle {input} -sid2 {wildcards.sample} -auto -outfile {output}"

rule viz_alignment:
    input:
        expand(
            path.join(
                "pairwise-alignments",
                "segment-{{segment}}-pairwise-alignments",
                "{sample}-segment-{{segment}}-alignment.needle"
            ),
            sample=SAMPLES
        )
    output:
        "visualizations/{segment}.pdf"
    script:
        "viz_alignment.py"

rule ref_segment_lengths:
    input:
        expand(
            path.join(
                "reference",
                "segment-{segment}-reference",
                "segment-{segment}-reference.fasta"
            ),
            segment=SEGMENTS
        )
    output:
        "reference-segment-lengths.json"
    script:
        "ref_segment_lengths.py"
