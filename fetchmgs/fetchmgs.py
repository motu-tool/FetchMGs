import argparse
from Bio import SeqIO
import glob
import re
import os
import subprocess
import sys
import pathlib
from shutil import which

def parse_cutoffs(args):
    """
    Parse the correct cutoffs file according to the run mode
    """
    cutoffs = {}
    with open(args.b, 'r') as fi:
        for line in fi:
            line = line.strip()
            if line[0] == '#':
                if args.mode == 'calibration':
                    if line != '#UNCALIBRATED CUTOFFS FILE':
                        sys.stderr.write(f'ERROR: When calibrating, the default bit score file provided with FetchMGs has to be used')
                        sys.exit(1)
                elif args.v:
                    if line != '#CALIBRATED CUTOFFS FILE - BEST HITS':
                        sys.stderr.write(f'ERROR: To extract using the option v|verybesthit, the file has to be calibrated with the -v option')
                        sys.exit(1)
                elif line != '#CALIBRATED CUTOFFS FILE - ALL HITS':
                    sys.stderr.write(f'ERROR: To extract without the option v|verybesthit, the file has to be calibrated without the -v option')
                    sys.exit(1)
            else:
               items = line.split('\t')
               cutoffs[items[0]] = int(items[1])
    return(cutoffs)

def run_hmmsearch(hmm, hmm_path, cutoff, args):
    '''
    Run hmmsearch from HMMER3, one hmm file against one set of protein records, then parse the results.
    '''
    hmmsearch_path = os.path.join(args.x, 'hmmsearch')
    out_path = os.path.join(args.o, f'hmmResults/{hmm}.out')
    tbl_path = os.path.join(args.o, f'hmmResults/{hmm}.dom')
    try:
        subprocess.run(f'hmmsearch --noali --notextw --cpu {args.t} -T {cutoff} -o {out_path} --domtblout {tbl_path} {hmm_path} {args.file}', shell=True, check=True)
    except subprocess.CalledProcessError as err:
        sys.stderr.write(f'\t{err}\n')
        sys.exit(1)
    hits = parse_hmmsearch(tbl_path)
    return(hits)

def parse_hmmsearch(tbl_path):
    '''
    Parse the standard output from hmmsearch to get scores by sequence.
    '''
    hits = {}
    with open(tbl_path, 'r') as fi:
        for line in fi:
            line = line.strip()
            if line[0] != '#':
                items = re.split('\s+', line)
                hits[items[0]] = float(items[7])
    return (hits)

def sort_hits(hits):
    '''
    Rearrange all hits from all hmms to find the best hmm for each sequence
    '''
    seqs = [seq for hmm in hits for seq in hits[hmm]]
    hits_by_seq = {seq:{hmm:hits[hmm][seq] for hmm in hits if seq in hits[hmm]} for seq in seqs}
    best_hits_by_seq = {seq:max(hits_by_seq[seq].items(), key=lambda x:x[1]) for seq in seqs}
    hits_by_hmm = {hmm:{seq:best_hits_by_seq[seq][1] for seq in seqs if hmm == best_hits_by_seq[seq][0]} for hmm in hits}

    return(seqs, hits_by_hmm)

def score_cutoff(pos, neg, nvalid, cutoff):
    '''
    Score a given cutoff based on how it would score known positives
    '''
    tp = len([x for x in pos if x >= cutoff])
    fp = len([x for x in neg if x >= cutoff])
    fn = nvalid - tp

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    fscore = 2 / ((1 / precision) + (1 / recall))

    return ([cutoff, tp, fp, fn, precision, recall, fscore])

def cli():
    parser = argparse.ArgumentParser(
        description='FetchMGs extracts the 40 single copy universal marker genes (decribed in Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes in an easy and accurate manner.',
        formatter_class=argparse.RawTextHelpFormatter)

    # Add -m option for legacy - it does nothing
    parser.add_argument('-m', '-mode', action='store_true', help='FetchMGs mode, see below')

    subparsers = parser.add_subparsers(title='modes', description='valid modes', dest='mode')
    subparsers.required = True
    PACKAGE_DIR = str(pathlib.Path(__file__).parent)
    ext_parser = subparsers.add_parser('extraction', help='extract marker genes from sequences',
                                       formatter_class=argparse.RawTextHelpFormatter)
    ext_parser.add_argument('file',
                            help='multi-FASTA file with protein sequences from which universal single-copy marker genes should be extracted')
    ext_parser.add_argument('-c', '-cog_used', nargs='+', default='all',
                            help='orthologous group id(s) to be extracted; example: "COG0012"')
    ext_parser.add_argument('-o', '-outdir', default='output', help='output directory')
    ext_parser.add_argument('-b', '-bitscore', default=None, help='path to bitscore cutoff file')
    ext_parser.add_argument('-l', '-library', default=os.path.join(PACKAGE_DIR, 'data'),
                            help='path to directory that contains hmm models')
    ext_parser.add_argument('-p', '-protein', action='store_true',
                            help='set if nucleotide sequences file for <protein sequences> is not available')
    ext_parser.add_argument('-d', '-dnaFastaFile',
                            help='multi-FASTA file with nucleotide sequences; not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes')
    ext_parser.add_argument('-v', '-verybesthit_only', action='store_true',
                            help='only extract the best hit of each COG from each genome\nrecommended to use, if extracting sequences from multiple reference genomes in the same file do not use it for metagenomes\nif this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have "project_id=XXX" in the header\nalternatively, set -i to ignore the headers; then, the best hit of each OG in the whole input file will be selected, regardless of the headers used')
    ext_parser.add_argument('-i', '-ignore_headers', action='store_true',
                            help='if this option is set in addition to -v, the best hit of each COG will be selected\nrecommended to use, if extracting sequences from a single genome in the same file')
    ext_parser.add_argument('-t', '-threads', default=1, help='number of processors/threads to be used')
    ext_parser.add_argument('-x', '-executable', default="",
                            help='path to executables used by this script\nif set to \'\', will search for executables in $PATH (default)')

    cal_parser = subparsers.add_parser('calibration',
                                       help='calibrate bitscores using results from extraction and a mapping file of known OGs',
                                       formatter_class=argparse.RawTextHelpFormatter)
    cal_parser.add_argument('file', help='file with sequences that include marker genes (true positives)')
    cal_parser.add_argument('map', help='tab-delimited file with true positive protein identifiers and COG ID')
    cal_parser.add_argument('-c', '-cog_used', nargs='+', default='all',
                            help='orthologous group id(s) to be extracted; example: "COG0012"')
    cal_parser.add_argument('-o', '-outdir', default='output', help='output directory')
    cal_parser.add_argument('-l', '-library', default=os.path.join(PACKAGE_DIR, 'data'),
                            help='path to directory that contains hmm models')
    cal_parser.add_argument('-p', '-protein', action='store_true',
                            help='set if nucleotide sequences file for <protein sequences> is not available')
    cal_parser.add_argument('-d', '-dnaFastaFile',
                            help='multi-FASTA file with nucleotide sequences; not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes')
    cal_parser.add_argument('-v', '-verybesthit_only', action='store_true',
                            help='only extract the best hit of each COG from each genome\nrecommended to use, if extracting sequences from multiple reference genomes in the same file do not use it for metagenomes\nif this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have "project_id=XXX" in the header\nalternatively, set -i to ignore the headers; then, the best hit of each OG in the whole input file will be selected, regardless of the headers used')
    cal_parser.add_argument('-i', '-ignore_headers', action='store_true',
                            help='if this option is set in addition to -v, the best hit of each COG will be selected\nrecommended to use, if extracting sequences from a single genome in the same file')
    cal_parser.add_argument('-t', '-threads', default=1, help='number of processors/threads to be used')
    cal_parser.add_argument('-x', '-executable', default="",
                            help='path to executables used by this script\nif set to \'\', will search for executables in $PATH (default)')

    args = parser.parse_args()
    return (args)

def import_files(args):
    # Compile list of HMM models
    hmms = glob.glob(f'{args.l}/*.hmm')
    if args.c == 'all':
        args.c = [os.path.splitext(os.path.split(x)[1])[0] for x in hmms]
    args.c = sorted(args.c)
    hmms = {x: f'{args.l}/{x}.hmm' for x in args.c}

    # Parse cutoff file
    if args.mode == 'calibration':
        args.b = f'{args.l}/MG_BitScoreCutoffs.uncalibrated.txt'

    if args.b is None:
        if args.v:
            args.b = f'{args.l}/MG_BitScoreCutoffs.verybesthit.txt'
        else:
            args.b = f'{args.l}/MG_BitScoreCutoffs.allhits.txt'
    cutoffs = parse_cutoffs(args)

    # Parse validated genes map if needed
    if args.mode == 'calibration':
        valid_map = {}
        for hmm in hmms.keys():
            valid_map[hmm] = []
        with open(args.map, 'r') as fi:
            for line in fi:
                line = line.strip()
                items = line.split("\t")
                valid_map[items[1]].append(items[0])
    else:
        valid_map = None

    # Set up protein sequence generators
    prot_records = SeqIO.parse(args.file, 'fasta')
    sys.stdout.write(f'Protein sequences: {args.file}\n')

    # Check for nucleotide file and set up generators
    if not args.p:
        if args.d is not None:
            nucl_records = SeqIO.parse(args.d, 'fasta')
            sys.stdout.write(f'Nucleotide sequences: {args.d}\n')
        else:
            nucl_file = f'{os.path.splitext(args.file)[0]}.fna'
            if os.path.exists(nucl_file):
                nucl_records = SeqIO.parse(nucl_file, 'fasta')
                sys.stdout.write(f'Nucleotide sequences: {nucl_file}\n')
            else:
                nucl_records = None
                sys.stdout.write(f'Nucleotide sequences: none found\n')
    else:
        nucl_records = None
        sys.stdout.write(f'Nucleotide sequences: none specified\n')

    return (hmms, cutoffs, valid_map, prot_records, nucl_records)

def extraction(args, hmms, cutoffs):
    # Go through HMMs and save results
    results = {}
    for hmm in hmms.keys():
        sys.stdout.write(f'    {hmm}\n')
        results[hmm] = run_hmmsearch(hmm, hmms[hmm], cutoffs[hmm], args)
    hit_ids, results = sort_hits(results)

    # Get tax_ids if there are multiple genomes
    if not args.i:
        if len(hit_ids[0].split(".")) > 1:
            tax_ids = set(x.split(".")[0] for x in hit_ids)
        else:
            tax_ids = None
    else:
        tax_ids = None 

    # Filter if -v
    if args.v:
        for hmm in hmms.keys():
            if tax_ids is not None:
                for tax_id in tax_ids:
                    genome_results = {k: v for k, v in results[hmm].items() if k.split(".")[0] == tax_id}
                    best_hit = max(genome_results, key=genome_results.get)
                    [results[hmm].pop(k) for k in genome_results.keys() if k != best_hit]
            else:
                best_hit = max(results[hmm], key=results[hmm].get)
                results[hmm] = {best_hit: results[hmm][best_hit]}

    return (results, hit_ids)

def calibration(args, results, valid_map, hmms, cutoffs):
    # Calibrate
    new_cutoffs = {}
    for hmm in hmms.keys():
        pos = sorted(v for k, v in results[hmm].items() if k in valid_map[hmm])
        neg = sorted(v for k, v in results[hmm].items() if k not in valid_map[hmm])
        max_cutoff = max(neg + pos)
        max_cutoff = int(max_cutoff/10)*10

        cutoff_scores = []
        min_cutoff = int(cutoffs[hmm]/10)*10
        for cutoff in range(min_cutoff, max_cutoff+10, 10):
            cutoff_scores.append(score_cutoff(pos, neg, len(valid_map[hmm]), cutoff))
        cutoff_scores = sorted(cutoff_scores, key=lambda x: (-x[6], -x[0]))
        new_cutoffs[hmm] = cutoff_scores[0]
    return (new_cutoffs)


def output_cutoffs(args, new_cutoffs, hmms):
    # Output bitscore file
    if args.v:
        header = '#CALIBRATED CUTOFFS FILE - BEST HITS\n'
        file = 'MG_BitScoreCutoffs.verybesthit.txt'
    else:
        header = '#CALIBRATED CUTOFFS FILE - ALL HITS\n'
        file = 'MG_BitScoreCutoffs.allhits.txt'
    with open(os.path.join(args.o, file), 'w') as fo:
        fo.write("#CALIBRATED CUTOFFS FILE - BEST HITS\n")
        for hmm in hmms.keys():
            fo.write(f'{hmm}\t{new_cutoffs[hmm][0]}\t{new_cutoffs[hmm][4]}\t{new_cutoffs[hmm][5]}\n')


def output_results(args, hmms, results, hit_ids, prot_records, nucl_records):
    # Only keep sequence records with hits
    prot_hits = {x.id: x for x in prot_records if x.id in hit_ids}
    if nucl_records is not None:
        nucl_hits = {x.id: x for x in nucl_records if x.id in hit_ids}

    # Output results table
    with open(os.path.join(args.o, "marker_genes_scores.tsv"), 'w') as fo:
        for hmm in hmms.keys():
            for id, score in sorted(results[hmm].items()):
                fo.write(f'{id}\t{score}\t{hmm}\t')
                project_id = re.search('project_id="(.+?)"', prot_hits[id].description)
                if project_id is not None:
                    tax_id = prot_hits[id].id.split(".")[0]
                    fo.write(f'{tax_id}.{project_id.group(1)}\n')
                else:
                    fo.write('-\n')

    # Output sequences
    for hmm in hmms.keys():
        with open(os.path.join(args.o, f'{hmm}.faa'), 'w') as fo:
            ids = sorted(results[hmm], key=results[hmm].get, reverse=True)
            for id in ids:
                seq = prot_hits[id]
                fo.write(f'>{seq.description}\n{seq.seq}\n')
        if nucl_records is not None:
            with open(os.path.join(args.o, f'{hmm}.fna'), 'w') as fo:
                ids = sorted(results[hmm], key=results[hmm].get, reverse=True)
                for id in ids:
                    seq = nucl_hits[id]
                    fo.write(f'>{seq.description}\n{seq.seq}\n')

def main():
    # Possible structural changes:
    #   Move makedirs to import_files (change name to setup or something)
    #   Move import of valid map to calibrate function

    args = cli()
    print(args)
    # Check that hmmsearch exists
    hmmsearch_path = os.path.join(args.x, 'hmmsearch')
    if which(hmmsearch_path) is None:
        sys.stderr.write(f'ERROR: hmmsearch cannot be found at {hmmsearch_path}.')
        sys.exit(1)

    # Create output directories before running HMMER
    os.makedirs(args.o, exist_ok = True)
    os.makedirs(os.path.join(args.o, 'hmmResults'), exist_ok = True)

    hmms, cutoffs, valid_map, prot_records, nucl_records = import_files(args)

    results, hit_ids = extraction(args, hmms, cutoffs)

    if args.mode == 'calibration':
        new_cutoffs = calibration(args, results, valid_map, hmms, cutoffs)
        output_cutoffs(args, new_cutoffs, hmms)

    output_results(args, hmms, results, hit_ids, prot_records, nucl_records)



if __name__ == '__main__':
    main()
