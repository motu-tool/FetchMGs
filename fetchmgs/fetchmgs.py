import argparse
import sys
import pathlib
import logging
import psutil
import pyrodigal
import pyhmmer
import Bio.SeqIO.FastaIO as FastaIO
import gzip
import collections
import tqdm
from typing import List

__author__ = ('Hans-Joachim Ruscheweyh (hansr@ethz.ch), '
              'Chris Field, '
              'Shinichi Sunagawa')
__version__ = '2.0.1'
__date__ = '10 Apr 2025'
__license__ = "GPL - v3"
__maintainer__ = "Hans-Joachim Ruscheweyh"
__email__ = 'hansr@ethz.ch'


def load_fetchmgs_files(very_best):
    data_folder = pathlib.Path(__file__).parent.absolute().joinpath('data')
    cutoffs_file = data_folder.joinpath('MG_BitScoreCutoffs.allhits.txt')
    cutoffs = {}
    if very_best:
        cutoffs_file = data_folder.joinpath('MG_BitScoreCutoffs.verybesthit.txt')

    with open(cutoffs_file, 'r') as f:
        f.readline()
        for line in f:
            splits = line.strip().split('\t')
            cutoffs[splits[0]] = float(splits[1])

    hmm_files = data_folder.glob('*hmm')
    cog_2_cutoff_hmm_file = {}
    for hmm_file in hmm_files:
        cog = str(hmm_file).split('/')[-1].split('.')[0]
        cutoff = cutoffs[cog]
        cog_2_cutoff_hmm_file[cog] = (cutoff, hmm_file)
    if len(cog_2_cutoff_hmm_file) != len(cutoffs):
        logging.error('There are marker genes in the cutoffs file for which hmm files are missing. Quitting ...')
        shutdown(1)
    return cog_2_cutoff_hmm_file


def get_file_handle(f):
    if is_gzipped(f):
        return gzip.open(f, 'rt')
    else:
        return open(f, 'r')


def check_sequence_type(input_files, expected_sequence_type):
    allowed_sequence_types = ['NT', 'AA']
    if expected_sequence_type not in allowed_sequence_types:
        logging.error(f'Expected sequence type {expected_sequence_type}. Quitting ...')
        shutdown(1)
    nucleotide_bases = set("ATGCUN")
    amino_acids = set("ARNDCQEGHILKMFPSTWYVX*")
    file_2_sequences = collections.defaultdict(list)

    for input_file in input_files:
        input_file_handle = get_file_handle(input_file)
        for (header, sequence) in FastaIO.SimpleFastaParser(input_file_handle):
            file_2_sequences[input_file].append(sequence)
            if len(file_2_sequences[input_file]) > 100:
                break
        input_file_handle.close()
    input_file_2_classification = {}
    input_file_2_average_length = {}
    for input_file, sequences in file_2_sequences.items():
        seqs = ''
        seq_lengths = []
        for seq in sequences:
            seq_clean = seq.upper().replace(" ", "")
            seqs = seqs + seq_clean
            seq_lengths.append(len(seq_clean))
        input_file_2_average_length[input_file] = sum(seq_lengths) / len(seq_lengths)
        chars = set(seq_clean)
        if chars.issubset(nucleotide_bases):
            input_file_2_classification[input_file] = "NT"
        elif chars.issubset(amino_acids):
            input_file_2_classification[input_file] = "AA"
        else:
            input_file_2_classification[input_file] = "UK"


    for input_file, sequence_type in input_file_2_classification.items():
        if sequence_type != expected_sequence_type:
            logging.error(f'Expected sequence type is {expected_sequence_type} but some sequences in {input_file} match sequence type {sequence_type}. Quitting ...')
            shutdown(1)

    return True






def extraction_genes(input_files: List[pathlib.Path], nucleotide_files: List[pathlib.Path], output_folder: pathlib.Path, threads: int, very_best: bool, genes_called: bool = False) -> None:
    """
    Extract the 40 marker genes from genes. Takes a list of protein
    files as input (input_files) and searches them against the HMM
    models where cutoffs have been calibrated.
    Writes the resulting marker genes in nucleotide and protein
    space to files together with the bitscores

    Args:
        input_files (List[pathlib.Path]): A list of input files containing genes in protein space. All files are required to exist.

        nucleotide_files (List[pathlib.Path]): A list of input files containing genes in nucleotide space. If no nucleotide files exist, then submit a list of None in the same length as the input files.

        output_folder pathlib.Path: A folder where the output files will be written. Will be created if it does not exist.

        threads (int): The number of threads to use for parallel processing.

        very_best (bool): By default all marker genes passing the cutoff will be reported. If this flag is set to True then only the best copy of the marker genes will be reported.

        genes_called (bool): False by default. Will be set to true if the extraction_genomes routine has also been executed

    Returns:
        None -> This method will write the results to output files and will not return any value
    """

    logging.info(f'Starting marker gene extraction from {len(input_files)} protein files.')
    output_folder.mkdir(exist_ok=True, parents=True)
    cog_2_cutoff_hmm_file = load_fetchmgs_files(very_best)
    hmms = []
    for cog, (cutoff, f) in cog_2_cutoff_hmm_file.items():
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
            hmm.cutoffs.trusted = (cutoff, cutoff)
            hmms.append(hmm)

    for prot_file, nucl_file in tqdm.tqdm(zip(input_files, nucleotide_files), total=len(input_files), unit=f'protein files'):

        if not pathlib.Path(prot_file).is_file():
            logging.error(f'Protein file {prot_file} does not exist. Quitting ...')
            shutdown(1)
        if nucl_file:
            if not pathlib.Path(nucl_file).is_file():
                logging.error(f'Nucleotide file {nucl_file} does not exist. Quitting ...')
                shutdown(1)

        proteins = None
        with pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seqs_file:
            proteins = seqs_file.read_block()

        # for each gene, collect the cog and score
        gene_2_cogs = collections.defaultdict(list)
        for hits in pyhmmer.hmmsearch(hmms, proteins, bit_cutoffs="trusted", cpus=threads):
            cog = hits.query.name.decode()
            for hit in hits:
                if hit.included:
                    gene_2_cogs[hit.name.decode()].append((cog, hit.score))

        # then pick for each gene the best scoring cog
        cog_2_hits = collections.defaultdict(list)
        for gene, cog_and_score in gene_2_cogs.items():
            (cog, score) = sorted(cog_and_score, key=lambda x: x[1], reverse=True)[0]
            cog_2_hits[cog].append((gene, score))

        # then sort by cog and see if very best hit is enabled
        if not very_best:
            cog_2_final_hits = cog_2_hits
        else:
            cog_2_final_hits = collections.defaultdict(list)
            for cog, hits in cog_2_hits.items():
                best_hit = sorted(hits, key=lambda x: x[1], reverse=True)[0]
                cog_2_final_hits[cog] = [best_hit]

        cog_2_final_hits.pop('COG0086', None)
        gene_2_cog_2_score = {}
        for cog, final_hits in cog_2_final_hits.items():
            for gene, score in final_hits:
                gene_2_cog_2_score[gene] = (cog, score)

        basename = str(prot_file).split('/')[-1]
        if genes_called:
            basename = basename.replace('.genes.faa', '')
            fna_out_file_name = output_folder.joinpath(basename + '.fetchMGs.fna')
            faa_out_file_name = output_folder.joinpath(basename + '.fetchMGs.faa')
            score_out_file_name = output_folder.joinpath(basename + '.fetchMGs.scores')
        else:
            fna_out_file_name = output_folder.joinpath(basename + '.fetchMGs.fna')
            faa_out_file_name = output_folder.joinpath(basename + '.fetchMGs.faa')
            score_out_file_name = output_folder.joinpath(basename + '.fetchMGs.scores')
        with open(faa_out_file_name, 'w') as faa_out_handle, open(score_out_file_name, 'w') as scores_out_handle:
            scores_out_handle.write('#protein_sequence_id\tHMM bit score\tCOG\n')
            seen_genes_nucl = set()
            seen_genes_prot = set()
            if nucl_file:
                with get_file_handle(nucl_file) as in_fna_handle, open(fna_out_file_name, 'w') as fna_out_handle:
                    for (header, sequence) in FastaIO.SimpleFastaParser(in_fna_handle):
                        header = header.split()[0]
                        if header in gene_2_cog_2_score:
                            cog, score = gene_2_cog_2_score[header]
                            fna_out_handle.write(f'>{header}.{cog}\n{sequence}\n')
                            seen_genes_nucl.add(header)
            with get_file_handle(prot_file) as in_faa_handle:
                for (header, sequence) in FastaIO.SimpleFastaParser(in_faa_handle):
                    header = header.split()[0]
                    if header in gene_2_cog_2_score:
                        cog, score = gene_2_cog_2_score[header]
                        faa_out_handle.write(f'>{header}.{cog}\n{sequence}\n')
                        scores_out_handle.write(f'{header}\t{int(score)}\t{cog}\n')
                        seen_genes_prot.add(header)
            if nucl_file:
                if seen_genes_nucl != seen_genes_prot:
                    logging.error('Some genes in the nucleotide file are missing. Quitting ...')
                    logging.error(seen_genes_prot.symmetric_difference(seen_genes_nucl))
                    shutdown(1)
    logging.info(f'Finished marker gene extraction.')


def extraction_genomes(input_files: List[pathlib.Path], output_folder: pathlib.Path, mode: str, threads: int, very_best: bool):
    """
    Extract genes from the input genomes/metagenomes and then also extract
    marker genes.

    Args:
        input_files (List[pathlib.Path]): A list of input files containing genomes or metagenomes. All files are required to exist.

        output_folder pathlib.Path: A folder where the output files will be written. Will be created if it does not exist.

        mode: str: can be either genome or metagenome. Decides on the mode in which the genes are being called.

        threads (int): The number of threads to use for parallel processing.

        very_best (bool): By default all marker genes passing the cutoff will be reported. If this flag is set to True then only the best copy of the marker genes will be reported.

    Returns:
        None -> This method will write the results to output files and will not return any value
    """
    logging.info(f'Starting gene calling from {len(input_files)} {mode} files.')
    # 1. extract the genes from the genomes/metagenomes
    # 2. run the extraction_genes routine
    output_folder.mkdir(exist_ok=True, parents=True)
    if mode == 'metagenome':
        gene_finder = pyrodigal.GeneFinder(meta=True, closed=True, mask=True)
    else:
        gene_finder = pyrodigal.GeneFinder(meta=False, closed=True, mask=True)

    nucleotide_files = []
    protein_files = []
    for input_file in tqdm.tqdm(input_files, total=len(input_files), unit=f'{mode}s', position=0):
        if not pathlib.Path(input_file).is_file():
            logging.error(f'{input_file} does not exist. Quitting ...')
            shutdown(1)

        scaffold_2_sequence = {}
        fname = str(input_file).split('/')[-1]
        fna_basename = f'{str(output_folder)}/{fname}.genes.fna'
        faa_basename = f'{str(output_folder)}/{fname}.genes.faa'
        nucleotide_files.append(pathlib.Path(fna_basename))
        protein_files.append(pathlib.Path(faa_basename))
        infile = None
        if str(input_file).endswith('.gz'):
            infile = gzip.open(input_file, 'rt')
        else:
            infile = open(input_file, 'r')
        for (header, sequence) in FastaIO.SimpleFastaParser(infile):
            header = header.split()[0]
            scaffold_2_sequence[header] = sequence
        if mode == 'genome':
            training_info = gene_finder.train(*(seq for seq in scaffold_2_sequence.values()))
        with open(fna_basename, 'w') as fi, open(faa_basename, 'w') as fo:
            for header, sequence in tqdm.tqdm(scaffold_2_sequence.items(), total=len(scaffold_2_sequence), unit='contigs', position=1, leave=False):
                genes = gene_finder.find_genes(sequence)
                genes.write_genes(fi, sequence_id=header)
                genes.write_translations(fo, sequence_id=header)
        infile.close()
    logging.info(f'Finished gene calling.')
    extraction_genes(protein_files, nucleotide_files, output_folder, threads, very_best, genes_called=True)

def is_gzipped(filename):
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'

def load_input_files_from_file(input_file):
    '''Checks input files and if they all exist.
    If the input file is a mapping file, return the contents

    '''

    input_handle = get_file_handle(input_file)
    is_seq_file = False
    for line in input_handle:
        if line.startswith('>'):
            is_seq_file = True
        break
    input_handle.close()
    tmp_input_files = []
    if is_seq_file:
        tmp_input_files = [input_file]
    else:

        broken_file = None
        input_handle = get_file_handle(input_file)
        for line in input_handle:
            if pathlib.Path(line.strip()).exists():
                tmp_input_files.append(pathlib.Path(line.strip()))
            else:
                broken_file = line.strip()
                break
        input_handle.close()

        if broken_file:
            logging.error(f'Input file is a mapping file. But some lines have non existing files. E.g. {broken_file}')
            shutdown(1)

    if len(tmp_input_files) == 0:
        logging.error(f'No valid input files found. Quitting ...')
        shutdown(1)
    return tmp_input_files

def parse_extraction():
    parser = argparse.ArgumentParser(usage=f'''Program: FetchMGs extracts the 40 
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

       -t INT         Number of threads. Default=[1]
       -v             Report only the very best hit per COG and input file. Only useful 
                        if input files contain genes from genomes or are genomes.
          ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    # Positional Parameters

    parser.add_argument(type=str, dest='input_file')
    parser.add_argument(type=str, choices=['gene', 'genome', 'metagenome'], dest='mode')
    parser.add_argument(type=str, dest='output_folder')

    # Input options
    parser.add_argument("-d", type=str, default=None, required=False)

    # Algorithm options

    parser.add_argument("-t", type=int, default=1)
    parser.add_argument("-v", action='store_true')

    startup()
    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        shutdown(1)

    # positionals

    input_file = pathlib.Path(args.input_file)
    output_folder = pathlib.Path(args.output_folder)
    mode = args.mode

    # optionals
    threads = args.t
    very_best = args.v
    nucleotide_file = args.d


    #checks
    input_sequence_files = load_input_files_from_file(input_file) # can be AA or NT files at this moment
    #nucleotide_input_files = [None] * len(input_sequence_files) # can be empty at this point

    if mode == 'gene':
        # check that provided sequences are AA
        check_sequence_type(input_sequence_files, 'AA')

        if args.d:
            nucleotide_input_files = load_input_files_from_file(nucleotide_file)
            check_sequence_type(nucleotide_input_files, 'NT')
            if len(nucleotide_input_files) != len(input_sequence_files):
                logging.error('Number of nucleotide files does not match number of protein files. Quitting ...')
                shutdown(1)
        else:
            nucleotide_input_files = [None] * len(input_sequence_files)



    elif mode == 'genome' or mode == 'metagenome':
        check_sequence_type(input_sequence_files, 'NT')
        if args.d:
            logging.error('-d parameter can\'t be set in metagenome or genome mode. Quitting ...')
            shutdown(1)
    else:
        logging.error(f'Invalid mode {mode}. Quitting ...')
        shutdown(1)


    if threads < 1:
        threads = 1
    if threads > psutil.cpu_count():
        logging.warning('Number of threads requested is above the number of CPUs.')

    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)



    if mode == 'gene':
        extraction_genes(input_sequence_files, nucleotide_input_files, output_folder, threads, very_best, genes_called=False)
    else:
        extraction_genomes(input_sequence_files, output_folder, mode, threads, very_best)


def shutdown(exitcode: int) -> None:
    """
    Securely shutdown fetchMGs.
    Args:
        exitcode: The exitcode

    Returns:
        None
    """
    #logging.info(f'fetchMGs shutting down with exitcode {exitcode}')
    sys.exit(exitcode)


def startup() -> None:
    """
    A method to group all functions that should be
    executed during startup of the fetchMGs tool.
    Returns:
        None
    """
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO, datefmt='%Y-%m-%d,%H:%M:%S')
    logging.info(f'fetchMGs {__version__} starting')


class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)



def main():
    parser = argparse.ArgumentParser(usage=f'''Program: FetchMGs extracts the 40 
    single copy universal marker genes (decribed in Ciccarelli et al., 
    Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes 
    in an easy and accurate manner.

    fetchMGs <command> [options]

      extraction     extract marker genes from sequences

    Type fetchMGs <command> to print the help menu for a specific command
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('command',
                        choices=["extraction"])
    args: argparse.Namespace = parser.parse_args(sys.argv[1:2])

    if args.command == 'extraction':
        parse_extraction()
    else:
        parser.print_usage()
        print(f'Unrecognized command {args}')
        shutdown(1)
    shutdown(0)


if __name__ == '__main__':
    main()

