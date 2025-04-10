# fetchMGs 2.0 for Python 3

fetchMGs (1.0 - 1.2) is copyright (c) 2019 Shinichi Sunagawa and Daniel R Mende.

fetchMGs (1.3) - written by Chris Field

fetchMGs (>=2.0) - written by Hans-Joachim Ruscheweyh


## Introduction
 
Phylogenetic markers are genes (and proteins) which can be used to reconstruct the phylogenetic history of different organisms. One classical phylogenetic marker is the 16S ribosomal RNA gene, which is often-used but is also known to be a sub-optimal phylogenetic marker for some organisms. Efforts to find a good set of protein coding phylogenetic marker genes (Ciccarelli et al., Science, 2006; Sorek et al., Science, 2007) lead to the identification of 40 universal single copy marker genes (MGs). These 40 marker genes occur in single copy in the vast majority of known organisms and they were used to successfully reconstruct a three domain phylogenetic tree (Ciccarelli et al., Science, 2006).

## What the software does
 
The program `fetchMGs` was written to extract the 40 MGs from genomes and metagenomes in an easy and accurate manner. This is done by utilizing Hidden Markov Models (HMMs) trained on protein alignments of known members of the 40 MGs as well as calibrated cutoffs for each of the 40 MGs. Please note that these cutoffs are only accurate when using complete protein sequences as input files. The output of the program are the protein sequences of the identified proteins, as well as their nucleotide sequences, if the nucleotide sequences of all complete genes are given as an additional input.


## Installation

FetchMGs and all its dependencies can be installed via `pip` and have been tested with Python 3.12.

```
$pip install fetchMGs
```



## Input

Users can submit genes in protein space or (from v2.0 on) longer nucleotide sequences from assembled genomes/metagenomes.

## Output

Per input sample (`SAMPLE`), fetchMGs will produce 3 output file:

1. `SAMPLE.fetchMGs.faa` --> the marker genes in protein space
2. `SAMPLE.fetchMGs.fna` --> the marker genes in nucleotide space
3. `SAMPLE.fetchMGs.scores` --> A link between marker genes and their bitscores

## Full program help


```
$fetchMGs

Program: FetchMGs extracts the 40
    single copy universal marker genes (decribed in Ciccarelli et al.,
    Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes
    in an easy and accurate manner.

    fetchMGs <command> [options]

      extraction     extract marker genes from sequences

    Type fetchMGs <command> to print the help menu for a specific command


```

### Extraction

```
$fetchMGs extraction

Program: FetchMGs extracts the 40
    single copy universal marker genes (decribed in Ciccarelli et al.,
    Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes
    in an easy and accurate manner.

    fetchMGs extraction [options]

    Positional arguments:
         FILE[ FILE]  Input file(s) - plain or gzipped. Can be either:
                            - 1-n genome assembly file(s), requires -m genome. Will
                                call genes before marker gene extraction.
                            - 1-n metagenome assembly file(s), requires -m metagenome. Will
                                call genes before marker gene extraction.
                            - 1-n gene file(s) in protein space, requires -m gene. nucleotide
                                sequences can be provided with -d parameter
                            - 1 text file with one line per input file. Requires 
                                -m parameter to enable "metagenome", "genome" or "gene" mode.
                                In "gene" mode another text file with samples in the
                                same order can be provided with -d parameter. 
    Input options:
       -d FILE[ FILE] Nucleotide files associated with protein files in -i. Same order as
                        files in -i required. Enabled only in -m gene mode. Can be either a 
                        list of files or a text file with one line per input file. 

    Output options:
       -o   FOLDER    Output folder for marker genes

    Algorithm options:
       -m STR         Mode of extraction Values: [gene, genome, metagenome]

       -t INT         Number of threads. Default=[1]
       -v             Report only the very best hit per COG and input file. Only useful
                        if input files contain genes from genomes or are genomes.

```

## Changelog

### 2.0.1

- Changed automatic detection of input files to `amino` for pyhmmer
- allow users to submit a file with a list of input files for positional and -d parameters


### 2.0.0

- Calibration mode was removed
- `hmmer` and `prodigal` were replaced with `pyhmmer` and `pyrodigal`
- Input is more flexible. Users can now submit multiple files and use different input formats:
	- Genes (`-m gene`)
	- Genomes (`-m genome`)
	- Metagenomes (`-m metagenome`)
- Output folder was cleaned up. Only one nucleotide and one protein file are generated compared to 40 in previous versions

### 1.3.0

- FetchMGs was ported from Perl to Python 3




