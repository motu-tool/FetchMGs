# FetchMGs 1.3 for Python 3

FetchMGs is copyright (c) 2019 Shinichi Sunagawa and Daniel R Mende
FetchMGs is released under the GNU General Public Licence v3.
Please see http://www.gnu.org/licenses/gpl.html and the seperately provided LICENSE file.

## Introduction
 
Phylogenetic markers are genes (and proteins) which can be used to reconstruct the phylogenetic history of different organisms. One classical phylogenetic marker is the 16S ribosomal RNA gene, which is often-used but is also known to be a sub-optimal phylogenetic marker for some organisms. Efforts to find a good set of protein coding phylogenetic marker genes (Ciccarelli et al., Science, 2006; Sorek et al., Science, 2007) lead to the identification of 40 universal single copy marker genes (MGs). These 40 marker genes occur in single copy in the vast majority of known organisms and they were used to successfully reconstruct a three domain phylogenetic tree (Ciccarelli et al., Science, 2006).

## What the software does
 
The program “FetchMGs” was written to extract the 40 MGs from genomes and metagenomes in an easy and accurate manner. This is done by utilizing Hidden Markov Models (HMMs) trained on protein alignments of known members of the 40 MGs as well as calibrated cutoffs for each of the 40 MGs. Please note that these cutoffs are only accurate when using complete protein sequences as input files. The output of the program are the protein sequences of the identified proteins, as well as their nucleotide sequences, if the nucleotide sequences of all complete genes are given as an additional input.

## Input

A fasta file with protein coding sequences, and optionally the gene sequences of the proteins. If the DNA sequences are available, the corresponding genes of the proteins, are also extracted.

## Output

The output of this software is saved within the specified output folder and consists of:
- 40 x COGxxxx.faa files (sequences of extracted proteins)
- 40 x COGxxxx.fna files (sequences of extracted genes)
- marker_genes_scores.table (tab-separated text with the columns protein, score, marker gene ID, genome identifier)
- hmmResults (specific output files from HMMer3)

## Full program help

usage: fetchmgs.py [-h] [-m] {extraction,calibration} ...

FetchMGs extracts the 40 single copy universal marker genes (decribed in Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes in an easy and accurate manner.

options:
  -h, --help            show this help message and exit
  -m, -mode             FetchMGs mode, see below

modes:
  valid modes

  {extraction,calibration}
    extraction          extract marker genes from sequences
    calibration         calibrate bitscores using results from extraction and a mapping file of known OGs

### Extraction

positional arguments:
  file                  multi-FASTA file with protein sequences from which universal single-copy marker genes should be extracted

options:
  -h, --help            show this help message and exit
  -c C [C ...], -cog_used C [C ...]
                        orthologous group id(s) to be extracted; example: "COG0012"
  -o O, -outdir O       output directory
  -b B, -bitscore B     path to bitscore cutoff file
  -l L, -library L      path to directory that contains hmm models
  -p, -protein          set if nucleotide sequences file for <protein sequences> is not available
  -d D, -dnaFastaFile D
                        multi-FASTA file with nucleotide sequences; not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes
  -v, -verybesthit_only
                        only extract the best hit of each COG from each genome
                        recommended to use, if extracting sequences from multiple reference genomes in the same file do not use it for metagenomes
                        if this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have "project_id=XXX" in the header
                        alternatively, set -i to ignore the headers; then, the best hit of each OG in the whole input file will be selected, regardless of the headers used
  -i, -ignore_headers   if this option is set in addition to -v, the best hit of each COG will be selected
                        recommended to use, if extracting sequences from a single genome in the same file
  -t T, -threads T      number of processors/threads to be used
  -x X, -executable X   path to executables used by this script
                        if set to '', will search for executables in $PATH (default)

### Calibration

positional arguments:
  file                  file with sequences that include marker genes (true positives)
  map                   tab-delimited file with true positive protein identifiers and COG ID

options:
  -h, --help            show this help message and exit
  -c C [C ...], -cog_used C [C ...]
                        orthologous group id(s) to be extracted; example: "COG0012"
  -o O, -outdir O       output directory
  -l L, -library L      path to directory that contains hmm models
  -p, -protein          set if nucleotide sequences file for <protein sequences> is not available
  -d D, -dnaFastaFile D
                        multi-FASTA file with nucleotide sequences; not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes
  -v, -verybesthit_only
                        only extract the best hit of each COG from each genome
                        recommended to use, if extracting sequences from multiple reference genomes in the same file do not use it for metagenomes
                        if this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have "project_id=XXX" in the header
                        alternatively, set -i to ignore the headers; then, the best hit of each OG in the whole input file will be selected, regardless of the headers used
  -i, -ignore_headers   if this option is set in addition to -v, the best hit of each COG will be selected
                        recommended to use, if extracting sequences from a single genome in the same file
  -t T, -threads T      number of processors/threads to be used
  -x X, -executable X   path to executables used by this script
                        if set to '', will search for executables in $PATH (default)

## Software dependencies

The fetchMGs script requires HMMER3 and Biopython:

HMMER3: http://www.hmmer.org/
Biopython: https://biopython.org/

