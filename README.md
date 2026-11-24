# fetchMGs 2 for Python 3

fetchMGs (1.0 - 1.2) is copyright (c) 2019 Shinichi Sunagawa and Daniel R Mende.

fetchMGs (1.3) - written by Chris Field

fetchMGs (>=2.0) - written by Hans-Joachim Ruscheweyh


## Introduction
 
Phylogenetic markers are genes (and proteins) which can be used to reconstruct the phylogenetic history of different organisms. One classical phylogenetic marker is the 16S ribosomal RNA gene, which is often-used but is also known to be a sub-optimal phylogenetic marker for some organisms. Efforts to find a good set of protein coding phylogenetic marker genes (Ciccarelli et al., Science, 2006; Sorek et al., Science, 2007) lead to the identification of 40 universal single copy marker genes (MGs). These 40 marker genes occur in single copy in the vast majority of known organisms and they were used to successfully reconstruct a three domain phylogenetic tree (Ciccarelli et al., Science, 2006).

## What the software does
 
The program `fetchMGs` was written to extract the 40 MGs from genomes and metagenomes in an easy and accurate manner. This is done by utilizing Hidden Markov Models (HMMs) trained on protein alignments of known members of the 40 MGs as well as calibrated cutoffs for each of the 40 MGs. Please note that these cutoffs are only accurate when using complete protein sequences as input files. The output of the program are the protein sequences of the identified proteins, as well as their nucleotide sequences, if the nucleotide sequences of all complete genes are given as an additional input.


## Installation

FetchMGs and all its dependencies can be installed via `pip` and have been tested with Python 3.10-3.14.

```
$pip install fetchMGs
```



## Input

Users can submit genes in protein space or (from v2.0 on) genomes/metagenomes in nucleotide space.

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

    fetchMGs extraction input_file mode output_folder[options]

    Positional arguments:
         FILE        Input file - plain or gzipped. Can be either:
                        - A genome assembly file (NT), requires -m genome. Will 
                            call genes before marker gene extraction.
                        - A metagenome assembly file (NT), requires -m metagenome. Will 
                            call genes before marker gene extraction.
                        - A gene file in protein space (AA), requires -m gene. nucleotide 
                            sequences can be provided with -d parameter
                        - A text file with one line per input file. Requires 
                            -m parameter to enable "metagenome", "genome" or "gene" mode.
                            In "gene" mode another text file with samples in the
                            same order can be provided with -d parameter. 
        
        STR          Mode of extraction Values: [gene, genome, metagenome]
        
        FOLDER       Output folder for marker genes
        
    Input options:
       -d FILE       Nucleotide file/Text file. Enabled only in the "gene" mode.
                     Requires same order of sequence files if submitted as
                     text file.
                     

    Algorithm options:

       -t INT         Number of threads, for gene calling and the marker gene
                        search. Default=[1]
       -b INT         Number of protein sequences searched at a time, across
                        input files. Memory grows with it. Default=[500000]
       -v             Report only the very best hit per COG and input file. Only useful 
                        if input files contain genes from genomes or are genomes.
          
```

## Changelog

### 2.2.0

- Same marker genes as 2.1.3; faster, and memory no longer grows with the input. On 6 threads:
  - 1M metagenome genes (gene mode): 201 s and up to 2.6 GB → 140 s and a flat 1.2 GB
  - 60 genomes (genome mode): 48 s → 24 s
  - a 153 Mbp metagenome (metagenome mode): 292 s → 90 s
- Marker gene search: every HMM is now searched by all threads, each on its share of the proteins (pyhmmer
  `parallel="targets"`). Before, each thread took whole HMMs, and every batch waited for COG0085, a quarter of the
  work, so no number of threads went past ~3.8x.
- Proteins are searched `-b` at a time (default 500,000), across input files: many small genomes share a batch,
  and every hit goes back to its file. pyhmmer's hits hold reference cycles that Python freed only rarely (GBs per
  million proteins); they are now collected after every batch.
- Gene calling uses the `-t` threads: genomes are called in parallel, each with its own Prodigal model, and the
  contigs of a metagenome are called in parallel and streamed, so the assembly is no longer held in memory.
  Metagenome mode calls every record; before, records with the same name were collapsed into one.
- The contig number in the `ID=` field of the `.genes.faa` and `.genes.fna` headers counts the contigs of each
  input file from 1, as Prodigal does and as a run with that file alone did. Before, it kept counting across the
  files of a run. The marker gene outputs carry only the gene names and are unchanged.
- Bug fixes:
  - Genomes and nucleotide gene files with IUPAC ambiguity codes (R, Y, K, M, S, W, B, D, H, V) were rejected as
    amino acids or unknown. They are accepted now; pyrodigal reads them as unknown bases, so called genes show N
    there and their proteins X.
  - Protein files with B, Z, J, U or O are accepted as amino acids.
  - An input file without sequences gives empty outputs and exit code 0, instead of "No valid input files found".
  - Results no longer depend on the file system or on the batches: the HMMs are loaded in sorted order (it decides
    the COG of a gene with two equal scores), and `-v` reports, of two genes with the same score, the one earlier
    in its file.
- Tests for threads, batches across files, ambiguous bases, empty input and ties.

### 2.1.3

- Updated `pyhmmer` to `~0.12.0`; fixed API compatibility (sequence name fields now return `str` instead of `bytes`)
- Loosened `pyrodigal` pin from `== 3.6.3.post1` to `>= 3.6.3.post1`
- Memory optimisation: protein sequences are now searched in chunks of 50,000 instead of loading entire files into memory at once
- Bug fixes:
  - Map files with empty or trailing blank lines no longer trigger false "file not found" errors
- Build system migrated to `pyproject.toml` only (removed `setup.py`, `setup.cfg`, `MANIFEST.in`)
- Regression test suite added (pytest + tox, GitHub Actions CI for Python 3.10–3.14)

### 2.1.0

- No functional updates = produces the same output
- Commandline interface refactored. 
  - fetchMGs expects 3 positional arguments, `input_file`, `mode` and `output_folder`
    - `input_file`: can be either a sequence file (genes in protein space, a genome or a metagenome) or a text file with one filepath per line to a genome/metagenome/protein file
    - `mode`: chose between `gene` (input file(s) are proteins), `genome` (input file(s) are genomes) or `metagenome` (input file(s) are metagenomes). 
    - `output_folder`: The output folder
- Features:
  - Checks if the input file is a map file or a sequence file now based on file content instead of file endings
  - Checks if the input (and nucleotide) file(s) are in the right format. E.g. a nucleotide file is represented by a different alphabet compared to amino acid sequences.

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





