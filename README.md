# H1N1 mutation count visualizer

Since the 1918 Spanish flu pandemic, there have been multiple H1N1 outbreaks,
but only two are commonly referred to as global pandemics:

* [1977 Russian flu pandemic][1977_russian_flu]
* [2009 swine flu pandemic][2009_swine_flu]

[1977_russian_flu]: https://en.wikipedia.org/wiki/1977_Russian_flu
[2009_swine_flu]: https://en.wikipedia.org/wiki/2009_swine_flu_pandemic

The objective of this pipeline is to provide a visual summary of mutations
across various H1N1 genomes. We hypothesize that by summarizing the number of
mutations along the length of these genomes, we will be able to visually observe
a unique pattern of mutation counts that distinguishes the Russian and swine flu
genomes from non-pandemic H1N1 genomes. If such a pattern is easily identifiable
upon visual inspection, we could use this pipeline to assess the pandemic
potential of future H1N1 outbreaks.

## How it works

Almost all cases of human H1N1 infections over the 20th and 21st centuries
were caused by descendents of the Spanish flu virus
([Taubenberger and Morens, 2006][mother_of_all_pandemics]). Therefore, to
identify mutations in the genome of H1N1 strains that arose following the
Spanish flu pandemic, this pipeline compares them with a reconstruction of the
Spanish flu genome. In this pipeline, we will refer to these post-Spanish flu
H1N1 strains as sample strains, and we will refer to the reconstructed Spanish
flu genome as the reference strain.

[mother_of_all_pandemics]: https://doi.org/10.3201/eid1201.050979

The sample strains this pipeline compares with the reference genome are:

| `{Influenza type}/{Location}/{Strain number}/{Year isolated}` | Pandemic strain? |
| --------------------------------------------------------- | ---------------- |
| A/PuertoRico/8/1934 | No |
| A/FortMonmouth/1/1947 | No |
| A/Albany/4836/1950 | No |
| A/NewJersey/1976 | No |
| A/USSR/92/1977 | **Yes** |
| A/India/6263/1980 | No |
| A/Taiwan/1/1986 | No |
| A/Texas/36/91 | No |
| A/Beijing/262/1995 | No |
| A/NewYork/222/2003 | No |
| A/California/07/2009| **Yes** |
| A/Alabama/13/2015 | No |
| A/Alabama/01/2020 | No |

The genome of H1N1 viruses consists of eight RNA molecules, called segments. As
shown in Fig. 1, we compare the sample H1N1 genomes to the reference genome on a
segment-by-segment basis to identify mutations in the sample genome. This
results in eight comparisons per sample strain, and since there are 13 sample
strains, this results in 8*13 = 104 comparisons in total.

![Overview of comparisons][comparisons-svg]

[comparisons-svg]: assets/images/comparisons.svg

### What happens during a comparison?

In each comparison, the fasta files for the reference and sample segment
sequences undergo a pairwise alignment using the
[Needleman-Wunsch algorithm][needle], which outputs alignment files in the
[`pairwise-alignments/`][pairwise-aligns] directory.

[pairwise-aligns]: pairwise-alignments/

Then, we iterate over the alignments using the [`viz_alignment.py`][viz-align]
script to mark the nucleotide position (relative to the reference genome) where
a mutation in the sample genome has occurred. Refer to Fig. 2. for an annotated
snippet from an alignment file, demonstrating how the file is parsed.

[viz-align]: viz_alignment.py

![Pairwise alignment example][pairwise-align]

[pairwise-align]: assets/images/pairwise-alignment.png

#### Note on sample size:

There is only one available genomic sequence per unique reference and
sample segment (i.e., _n_ = 1 in each comparison). Many of the strains in
this pipeline only had one complete genomic sequence available. I retrieved the
strains from [fludb.org][fludb].

[fludb]: https://www.fludb.org/brc/home.spg?decorator=influenza

### What happens after all comparisons are completed?

After counting the mutations from the pairwise alignment outputs, the
[`viz_alignment.py`][viz-align] script proceeds to generate three heatmaps for
each segment, encoding the amount of SNP, insertion, and deletion mutations per
100 nucleotide positions found in each sample strain, along the length of the
segment. The pandemic strain labels are colored red for clarity. The heatmaps
generated for segment 4 are shown in Fig. 3. as an example.

![Segment 4 heatmaps][heatmaps]

[heatmaps]: assets/images/heatmaps.svg

[`assets/reference-segment-lengths.json`][ref-segments-len] is used when
determining the number of bins needed when grouping each segment genome into
collections of 100 nucleotide bases. This file will be automatically generated
by the pipeline, if it goes missing, or the reference fasta files are updated.

[ref-segments-len]: assets/reference-segment-lengths.json

## Reference genome data

The reference genome segment fasta files can be found nested in the
[`reference/`][ref] directory. Full accession details can be found in the header
of each fasta file, including genbank accession number.

[ref]: reference/

Segments 1-3 and 5-8 are taken from the A/Brevig Mission/1/1918 Spanish flu
strain, which was reconstructed in the late 1990s from the frozen remains of a
woman that succumbed to the Spanish flu in 1918. The last 480 nucleotides of
segment 4 are missing from this particular strain, so segment 4 is taken from
another partially reconstructed strain: A/South Carolina/1/18. This approach of
combining the strains is similar to one taken by
[Carter and Sanford (2012)][reconstruction], who observed that segment 4 of the
two strains was nearly identical prior to the last 480 nucleotides missing in
A/Brevig Mission/1/1918.

[reconstruction]: https://doi.org/10.1186/1742-4682-9-42

## Sample genome data

The individual non-reference sample segment fasta files can be found nested in
the [`samples/`][samples] directory. Full accession details can be found in the
header of each fasta file, including genbank accession number.

Sample fasta files must be placed in the appropriate nested segment directory,
and observe the following naming scheme:

[samples]: samples/

`{sample}-segment-{segment}.fasta`

For example, the eighth segment fasta file for A/Alabama/01/2020 is saved in the
[samples/segment-8-samples][segment-8] directory as
`A-Alabama-01-2020-segment-8.fasta`. We replaced the `/` characters with `-`
characters, because you cannot use `/` in file names. The sample name is
flexible, primarily effecting the label displayed downstream on the heatmap.

[segment-8]: samples/segment-8-samples

For each sample, you must have files for all eight segments in the appropriate
directories.

### Configuration files

There are two additional, customizable files that can help configure the
pipeline. [`assets/samples.json`][samples-json] includes a list of all samples
that will be compared with the reference genome in the pipeline, so you do not
have to generate comparisons for every single sample file in the [samples/]
[samples] directory--only the samples you list in
[`assets/samples.json`][samples-json].

[samples-json]: assets/samples.json

[`assets/pandemic-strains.json`][pandemic-strains] specifies the samples that
belong the 1977 and 2009 pandemics. The main purpose of this file is to color
the pandemic strain labels in the visualized heatmaps red.

[pandemic-strains]: assets/pandemic-strains.json

## Results

The results of the pipeline are in the [`visualizations/`][viz] directory. Each
segment has a pdf file with a SNP, insertion, and deletion mutation heatmap.

[viz]: visualizations/

It is difficult to visually observe any noticeable mutation count trends that
distinguish the two pandemic sample strains from the non-pandemic sample
strains, in any segment. The number of mutations does seem to have increased
over the years, with the highest number of mutations found in the most recent
strains, but this increase seems to occur at relatively constant rate across
pandemic and non-pandemic strains. Thus, the pipeline results do not immediately
support our hypothesis that a visual summary of mutation counts along the length
of the genome will enable us to observe a trend unique to pandemic strains.

However, it may be beneficial to continue modifying the pipeline for further
testing of our hypothesis. Heatmaps have the potential encode thousands of data
points, while still allowing users to distinguish visual patterns. We could
potentially increase the granularity of the heatmap, and stop binning the
results per every 100 nucleotides. We could also include many more strains--one
from every year since 1918, or perhaps multiple strains per year.
[`viz_alignment.py`][viz-align] would need to be adjusted slightly, and we would
need to download more sample data. We may also need to use a more powerful
machine to generate our results in a timely manner.

# Installation

The pipeline is run using snakemake, in a conda environment. The dependencies
required by the conda environment are listed in the
[`environment.yaml`][env] file. Here is a description of the key dependencies:

[env]: environment.yaml

| Dependency | Version | Function |
| ---------- | ------- | -------- |
| snakemake-minimal | 6.10.0 | Lightweight snakemake installation. Used to run pipeline |
| python | 3.10.0 | Required by snakemake, and used to run several scripts. We specify a specific version because the functionality of certain Python operations has changed slightly over the years. |
| emboss | 6.6.0 | Provides the needle package used to run the [Needleman-Wunsch algorithm][needle] for pairwise alignments. |
| biopython | 1.79 | Provides AlignIO and SeqIO packages that are used to parse outputted alignment files and fasta sequence files respectively. |
| matplotlib-base | 3.4.3 | Used to generate visualizations. |

[needle]: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

## Install conda

If you are running a linux or Mac OS, you can check if you already have conda
installed by typing the following in the terminal:

`$ conda --version`

If you have windows, you can search for the conda prompt in your program files.

If you do not have conda installed, you can follow the instructions for
installing Miniconda, which in my opinion is the easiest and most efficient way
of installing conda. [You can find installation instructions for different
operating systems here.][miniconda]

[miniconda]: https://bit.ly/3HTXiq0

## Install git (Windows users only)

If you are running a Windows OS, you may need to install git. [Follow the
instructions for installing git on Windows here.][git]

[git]: https://bit.ly/3xw8SD2

## Install the pipeline

**Note for windows users: although the following commands are geared towards
Linux and Mac users, you should be able to open the git bash program and run the
following commands.**

**Note for first time git users: If you're having problems cloning the
repository, you can also click the green Code button at the top of the repo to
download and extract the repository as a ZIP file. Then, just navigate into the
extracted directory with the `cd` command.**

Clone this repository, and navigate into the cloned directory.

`$ git clone git@github.com:ivansg44/biof-501-project.git`

`$ cd biof-501-project`

Create the conda environment.

`$ conda env create --file environment.yaml`

And that's it! The pipeline should be ready to go.

# Usage

To run the pipeline, you must navigate to the cloned directory and activate the
conda environment.

`$ cd biof-501-project`

`$ conda activate biof-501-project`

Then, you can simply run:

`$ snakemake --cores`

And that's it. This command will use all your computer cores to run the
pipeline, if necessary.

***Warning: there are 104 pairwise alignments. This pipeline takes ~50 mins to
run on my machine, which has 4 cores. For the sake of time, you can remove
samples from [`assets/samples.json`][samples-json], so you only run the pipeline
on a subset of samples. If you need to cancel an already running workflow, hit
Ctrl+C.**

# Author

[@ivansg44][github] - This is my term project for BIOF 501.

[github]: https://github.com/ivansg44

# License

[MIT](LICENSE)
