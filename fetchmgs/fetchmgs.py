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
import concurrent.futures
import gc
import io
import re
import tqdm
from typing import List
from importlib.metadata import version

__author__ = ('Hans-Joachim Ruscheweyh (hansr@ethz.ch), '
              'Chris Field, '
              'Shinichi Sunagawa')
__version__ = version("fetchMGs")
__date__ = '23 Nov 2026'

# Number of protein sequences searched against the HMMs at a time (-b). The proteins of several files can share a
# batch, and a large file is split over several. Memory grows with the batch, not with the input.
DEFAULT_BATCH_SIZE = 500000

# Sequence type check (first 101 sequences of a file). A nucleotide file may hold IUPAC ambiguity codes; pyrodigal
# reads every letter but A, C, G and T as an unknown base, written as N in the genes. At least 90% of the letters must
# be A, C, G, T, U or N, so that a few short peptides made of nucleotide letters are not taken for DNA.
NUCLEOTIDE_CODES = set('ACGTUNRYKMSWBDHV')
PLAIN_NUCLEOTIDES = 'ACGTUN'
MIN_PLAIN_NUCLEOTIDE_SHARE = 0.9
AMINO_ACIDS = set('ARNDCQEGHILKMFPSTWYVX*BZJUO')

# At most this many contigs per thread wait for gene calling or writing in metagenome mode.
CONTIGS_PER_THREAD = 4

# The contig number in the ID=<contig>_<gene> field of pyrodigal's gene headers.
GENE_ID = re.compile(r' # ID=\d+_')
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

    # sorted: the HMM order decides which COG a gene gets when two score the same, so it must not depend on the
    # file system's listing order
    hmm_files = sorted(data_folder.glob('*hmm'))
    cog_2_cutoff_hmm_file = {}
    for hmm_file in hmm_files:
        cog = hmm_file.stem
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


def classify_sequences(sequences):
    """
    Classify sequences as nucleotides (NT), amino acids (AA) or unknown (UK).

    NT: only IUPAC nucleotide codes, at least MIN_PLAIN_NUCLEOTIDE_SHARE of them A, C, G, T, U or N.
    AA: only amino acid letters, including B, Z, J, U, O, X and the stop *.
    Checked in this order, as every nucleotide code but U is also an amino acid letter.

    Returns:
        'NT', 'AA', 'UK', or None if the sequences hold no letters at all
    """
    chars = set()
    plain = 0
    total = 0
    for sequence in sequences:
        sequence = sequence.upper().replace(" ", "")
        chars.update(sequence)
        plain += sum(sequence.count(base) for base in PLAIN_NUCLEOTIDES)
        total += len(sequence)
    if total == 0:
        return None
    if chars <= NUCLEOTIDE_CODES and plain >= MIN_PLAIN_NUCLEOTIDE_SHARE * total:
        return "NT"
    if chars <= AMINO_ACIDS:
        return "AA"
    return "UK"


def check_sequence_type(input_files, expected_sequence_type):
    allowed_sequence_types = ['NT', 'AA']
    if expected_sequence_type not in allowed_sequence_types:
        logging.error(f'Expected sequence type {expected_sequence_type}. Quitting ...')
        shutdown(1)

    for input_file in input_files:
        sequences = []
        with get_file_handle(input_file) as input_file_handle:
            for (_, sequence) in FastaIO.SimpleFastaParser(input_file_handle):
                sequences.append(sequence)
                if len(sequences) > 100:
                    break
        sequence_type = classify_sequences(sequences)
        if sequence_type is not None and sequence_type != expected_sequence_type:
            logging.error(f'Expected sequence type is {expected_sequence_type} but some sequences in {input_file} match sequence type {sequence_type}. Quitting ...')
            shutdown(1)

    return True






def load_hmms(very_best):
    hmms = []
    for cog, (cutoff, f) in load_fetchmgs_files(very_best).items():
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
            hmm.cutoffs.trusted = (cutoff, cutoff)
            hmms.append(hmm)
    return hmms


def search_marker_genes(protein_files, hmms, threads, chunk_size):
    """
    Search the proteins of all files against the HMMs, chunk_size proteins at a time. A batch may hold the proteins
    of several files (small genomes fill one together), and a large file is split over several batches. For the
    search each protein is named <file index>:<position in its file>:<name>, so genes of the same name in different
    files stay apart and every hit goes back to its file.

    Args:
        protein_files (List[pathlib.Path]): Files with genes in protein space. An empty file has no proteins.

        hmms (list): The HMMs, with their cutoffs as trusted cutoffs.

        threads (int): The number of threads to use for the search.

        chunk_size (int): The number of proteins searched at a time.

    Returns:
        Per file a dict gene -> [(cog, score), ...] in HMM order, and per file a dict gene -> position in its file
        (for the genes with a hit only)
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    file_2_gene_2_cogs = [collections.defaultdict(list) for _ in protein_files]
    file_2_gene_2_position = [{} for _ in protein_files]

    def search(batch):
        block = pyhmmer.easel.DigitalSequenceBlock(alphabet, batch)
        # parallel="targets": every HMM is searched by all threads, each on its share of the batch (balanced by
        # residues). By default each thread takes whole HMMs, and COG0085 alone is a quarter of the work.
        for hits in pyhmmer.hmmsearch(hmms, block, bit_cutoffs="trusted", cpus=threads, parallel="targets"):
            cog = hits.query.name
            for hit in hits:
                if hit.included:
                    file_index, position, gene = hit.name.split(':', 2)
                    file_index = int(file_index)
                    file_2_gene_2_cogs[file_index][gene].append((cog, hit.score))
                    file_2_gene_2_position[file_index][gene] = int(position)

    batch = []
    for file_index, prot_file in enumerate(tqdm.tqdm(protein_files, total=len(protein_files), unit='protein files')):
        try:
            seqs_file = pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=alphabet, format='fasta')
        except EOFError:  # an empty file: no proteins
            continue
        with seqs_file:
            for position, seq in enumerate(seqs_file):
                seq.name = f'{file_index}:{position}:{seq.name}'
                batch.append(seq)
                if len(batch) == chunk_size:
                    search(batch)
                    batch = []
                    # pyhmmer's hits hold reference cycles around large C buffers. Only a full garbage collection
                    # frees them, and Python runs one rarely: without this, memory grew by GBs per million proteins.
                    gc.collect()
    if batch:
        search(batch)
        batch = []
        gc.collect()
    return file_2_gene_2_cogs, file_2_gene_2_position


def extraction_genes(input_files: List[pathlib.Path], nucleotide_files: List[pathlib.Path], output_folder: pathlib.Path, threads: int, very_best: bool, genes_called: bool = False, chunk_size: int = DEFAULT_BATCH_SIZE) -> None:
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

        chunk_size (int): The number of protein sequences searched at a time, across files.

    Returns:
        None -> This method will write the results to output files and will not return any value
    """

    logging.info(f'Starting marker gene extraction from {len(input_files)} protein files.')
    output_folder.mkdir(exist_ok=True, parents=True)
    hmms = load_hmms(very_best)

    for prot_file, nucl_file in zip(input_files, nucleotide_files):
        if not pathlib.Path(prot_file).is_file():
            logging.error(f'Protein file {prot_file} does not exist. Quitting ...')
            shutdown(1)
        if nucl_file:
            if not pathlib.Path(nucl_file).is_file():
                logging.error(f'Nucleotide file {nucl_file} does not exist. Quitting ...')
                shutdown(1)

    # for each gene, collect the cog and score
    file_2_gene_2_cogs, file_2_gene_2_position = search_marker_genes(input_files, hmms, threads, chunk_size)

    for file_index, (prot_file, nucl_file) in enumerate(zip(input_files, nucleotide_files)):
        gene_2_cogs = file_2_gene_2_cogs[file_index]
        gene_2_position = file_2_gene_2_position[file_index]

        # then pick for each gene the best scoring cog (on equal scores the first in HMM order)
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
                # on equal scores the gene that comes first in the input file, whatever the batches were
                best_hit = sorted(hits, key=lambda x: (-x[1], gene_2_position[x[0]]))[0]
                cog_2_final_hits[cog] = [best_hit]

        cog_2_final_hits.pop('COG0086', None)
        gene_2_cog_2_score = {}
        for cog, final_hits in cog_2_final_hits.items():
            for gene, score in final_hits:
                gene_2_cog_2_score[gene] = (cog, score)

        basename = pathlib.Path(prot_file).name
        if genes_called:
            basename = basename.replace('.genes.faa', '')
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


def call_genes_genome(input_file, fna_file, faa_file):
    """
    Genome mode: train Prodigal on all contigs of the genome, then call the genes of each contig. Every genome gets
    its own GeneFinder, so genomes can be called in parallel threads (pyrodigal releases the GIL while it works).
    """
    scaffold_2_sequence = {}
    with get_file_handle(input_file) as infile:
        for (header, sequence) in FastaIO.SimpleFastaParser(infile):
            header = header.split()[0]
            scaffold_2_sequence[header] = sequence
    gene_finder = pyrodigal.GeneFinder(meta=False, closed=True, mask=True)
    if scaffold_2_sequence:
        gene_finder.train(*(seq for seq in scaffold_2_sequence.values()))
    with open(fna_file, 'w') as fi, open(faa_file, 'w') as fo:
        for header, sequence in scaffold_2_sequence.items():
            genes = gene_finder.find_genes(sequence)
            genes.write_genes(fi, sequence_id=header)
            genes.write_translations(fo, sequence_id=header)


def call_genes_metagenome(input_file, fna_file, faa_file, gene_finder, pool, threads):
    """
    Metagenome mode: the contigs are read one after the other and their genes called in parallel threads, all with
    the same meta-mode GeneFinder (it is never trained, so it can be shared). Genes are written in contig order, and
    at most CONTIGS_PER_THREAD contigs per thread are held in memory, not the whole assembly.
    """
    pending = collections.deque()
    with get_file_handle(input_file) as infile, open(fna_file, 'w') as fi, open(faa_file, 'w') as fo, \
            tqdm.tqdm(unit='contigs', position=1, leave=False) as progress:

        def write_next():
            contig_number, header, future = pending.popleft()
            genes = future.result()
            # The ID=<contig>_<gene> field: pyrodigal numbers the contigs of a shared GeneFinder in the order the
            # threads happen to start them, so it is set here to the contig's position in its file, as Prodigal
            # numbers the contigs of one input file. Contig names hold no whitespace, so ' # ID=' is only the field.
            for write, handle in ((genes.write_genes, fi), (genes.write_translations, fo)):
                buffer = io.StringIO()
                write(buffer, sequence_id=header)
                handle.write(GENE_ID.sub(f' # ID={contig_number}_', buffer.getvalue()))
            progress.update()

        for contig_number, (header, sequence) in enumerate(FastaIO.SimpleFastaParser(infile), start=1):
            pending.append((contig_number, header.split()[0], pool.submit(gene_finder.find_genes, sequence)))
            if len(pending) >= CONTIGS_PER_THREAD * threads:
                write_next()
        while pending:
            write_next()


def extraction_genomes(input_files: List[pathlib.Path], output_folder: pathlib.Path, mode: str, threads: int, very_best: bool, chunk_size: int = DEFAULT_BATCH_SIZE):
    """
    Extract genes from the input genomes/metagenomes and then also extract
    marker genes.

    Args:
        input_files (List[pathlib.Path]): A list of input files containing genomes or metagenomes. All files are required to exist.

        output_folder pathlib.Path: A folder where the output files will be written. Will be created if it does not exist.

        mode: str: can be either genome or metagenome. Decides on the mode in which the genes are being called.

        threads (int): The number of threads to use for parallel processing: genomes are called in parallel, the
            contigs of a metagenome are called in parallel, and the marker gene search uses them too.

        very_best (bool): By default all marker genes passing the cutoff will be reported. If this flag is set to True then only the best copy of the marker genes will be reported.

        chunk_size (int): The number of protein sequences searched at a time, across files.

    Returns:
        None -> This method will write the results to output files and will not return any value
    """
    logging.info(f'Starting gene calling from {len(input_files)} {mode} files.')
    # 1. extract the genes from the genomes/metagenomes
    # 2. run the extraction_genes routine
    output_folder.mkdir(exist_ok=True, parents=True)
    for input_file in input_files:
        if not pathlib.Path(input_file).is_file():
            logging.error(f'{input_file} does not exist. Quitting ...')
            shutdown(1)

    nucleotide_files = [output_folder / f'{pathlib.Path(input_file).name}.genes.fna' for input_file in input_files]
    protein_files = [output_folder / f'{pathlib.Path(input_file).name}.genes.faa' for input_file in input_files]
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        if mode == 'metagenome':
            gene_finder = pyrodigal.GeneFinder(meta=True, closed=True, mask=True)
            for input_file, fna_file, faa_file in tqdm.tqdm(zip(input_files, nucleotide_files, protein_files), total=len(input_files), unit=f'{mode}s', position=0):
                call_genes_metagenome(input_file, fna_file, faa_file, gene_finder, pool, threads)
        else:
            futures = [pool.submit(call_genes_genome, input_file, fna_file, faa_file)
                       for input_file, fna_file, faa_file in zip(input_files, nucleotide_files, protein_files)]
            for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures), unit=f'{mode}s', position=0):
                future.result()
    logging.info(f'Finished gene calling.')
    extraction_genes(protein_files, nucleotide_files, output_folder, threads, very_best, genes_called=True, chunk_size=chunk_size)

def is_gzipped(filename):
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'

def load_input_files_from_file(input_file):
    '''Checks input files and if they all exist.
    If the input file is a mapping file, return the contents

    '''

    first_line = None
    with get_file_handle(input_file) as input_handle:
        for line in input_handle:
            if line.strip():
                first_line = line
                break
    if first_line is None:
        # nothing but blank lines: a sequence file without sequences, which gives empty outputs
        logging.warning(f'{input_file} holds no sequences.')
        return [input_file]
    is_seq_file = first_line.startswith('>')
    tmp_input_files = []
    if is_seq_file:
        tmp_input_files = [input_file]
    else:

        broken_file = None
        with get_file_handle(input_file) as input_handle:
            for line in input_handle:
                if not line.strip():
                    continue
                if pathlib.Path(line.strip()).exists():
                    tmp_input_files.append(pathlib.Path(line.strip()))
                else:
                    broken_file = line.strip()
                    break

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
    
    Version: {__version__}

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
                        input files. Memory grows with it. Default=[{DEFAULT_BATCH_SIZE}]
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
    parser.add_argument("-b", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("-v", action='store_true')

    if not sys.argv[2:]:
        parser.print_usage()
        shutdown(0)

    startup()
    args = parser.parse_args(sys.argv[2:])

    # positionals

    input_file = pathlib.Path(args.input_file)
    output_folder = pathlib.Path(args.output_folder)
    mode = args.mode

    # optionals
    threads = args.t
    batch_size = args.b
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



    elif mode in ('genome', 'metagenome'):
        check_sequence_type(input_sequence_files, 'NT')
        if args.d:
            logging.error('-d parameter can\'t be set in metagenome or genome mode. Quitting ...')
            shutdown(1)


    if threads < 1:
        threads = 1
    if batch_size < 1:
        logging.error(f'-b must be at least 1 (got {batch_size}). Quitting ...')
        shutdown(1)
    if threads > psutil.cpu_count():
        logging.warning('Number of threads requested is above the number of CPUs.')

    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)



    if mode == 'gene':
        extraction_genes(input_sequence_files, nucleotide_input_files, output_folder, threads, very_best, genes_called=False, chunk_size=batch_size)
    else:
        extraction_genomes(input_sequence_files, output_folder, mode, threads, very_best, chunk_size=batch_size)


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
    
    Version: {__version__}
    
    fetchMGs <command> [options]

      extraction     extract marker genes from sequences

    Type fetchMGs <command> to print the help menu for a specific command
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('command',
                        choices=["extraction"])

    if not sys.argv[1:]:
        parser.print_usage()
        shutdown(0)

    args: argparse.Namespace = parser.parse_args(sys.argv[1:2])

    if args.command == 'extraction':
        parse_extraction()
    shutdown(0)


if __name__ == '__main__':
    main()

