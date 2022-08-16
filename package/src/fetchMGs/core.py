import argparse
from Bio import SeqIO
import glob
import re
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='fetchMGs extracts the 40 single copy universal marker genes (decribed in Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes in an easy and accurate manner.')

subparsers = parser.add_subparsers(title='modes', description='valid modes', dest='mode')
subparsers.required = True

ext_parser = subparsers.add_parser('extraction', help='extract marker genes from sequences')
ext_parser.add_argument('files', nargs="+", help='files with sequences from which universal single-copy marker genes should be extracted')
ext_parser.add_argument('-c', '--cogs', nargs='+', default='all', help='orthologous group id(s) to be extracted; example: "COG0012"')
ext_parser.add_argument('-o', '--outdir', default='fetchmgs', help='output directory')
ext_parser.add_argument('-b', '--bitscore', default='data/MG_BitScoreCutoffs.verybesthit.txt', help='path to bitscore cutoff file')
ext_parser.add_argument('-l', '--library', default='data', help='path to directory that contains hmm models')
ext_parser.add_argument('-a', '--all', action='store_true', help='output all results; default output is the best result for each COG in each file')
ext_parser.add_argument('-t', '--threads', default=1, help='number of processors/threads to be used')
ext_parser.add_argument('--hmmsearch', default='hmmsearch', help='path to hmmsearch executable if not in PATH')

cal_parser = subparsers.add_parser('calibration', help='calibrate bitscores using results from extraction and a mapping file of known OGs')
cal_parser.add_argument('file', help='file with sequences that include marker genes (true positives)')
cal_parser.add_argument('map', help='tab-delimited file with true positive protein identifiers and COG ID')
cal_parser.add_argument('-o', '--outdir', default='fetchmgs', help='output directory')
cal_parser.add_argument('-b', '--bitscore', help='path to bitscore cutoff file')
# None of this mode implemented yet

args = parser.parse_args()

# Amino acid codes that aren't ACGT
aas = ['P', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S']

# Compile list of HMM models
hmms = glob.glob(f'{args.library}/*.hmm')
if args.cogs == 'all':
    args.cogs = [os.path.splitext(os.path.split(x)[1])[0] for x in hmms]
hmms = {x:f'{args.library}/{x}.hmm' for x in args.cogs}

# Parse cutoff file
cutoffs = {}
with open(args.bitscore, 'r') as fi:
    for line in fi.readlines():
        line = line.strip()
        if line[0] != '#':
           items = line.split('\t')
           cutoffs[items[0]] = int(items[1])

def is_nucl(record):
    '''
    Guess whether a SeqRecord is protein or nucleotide sequence.
    '''
    if any(item in aas for item in record.seq.upper()):
        return(False)
    else:
        return(True)

def run_hmmsearch(hmm_path, stream, cutoff, threads, hmmsearch):
    '''
    Run hmmsearch from HMMER3, one hmm file against one set of protein records, then parse the results.
    '''
    try:
        output = subprocess.check_output(f'{hmmsearch} --noali --notextw --cpu {threads} -T {cutoff} {hmm_path} -', input=stream.encode(), shell=True).decode()
    except subprocess.CalledProcessError as err:
        sys.stdout.write(f'\t{err}\n')
        sys.exit(1)
    results = parse_hmmsearch(output)
    return(results)

def parse_hmmsearch(output):
    '''
    Parse the standard output from hmmsearch to get scores by sequence.
    '''
    hits = {}
    countdown = 1000
    for line in output.split("\n"):
        line = line.strip()
        if line[0:6] == 'Scores':
            countdown = 4
        if len(line) == 0:
            countdown = 1000
        if countdown == 0:
            items = re.split("\s+", line)
            hits[items[8]] = float(items[1])
        else:
            countdown -= 1
    return(hits)

# Main loop through files to create protein record stream to send to hmmsearch
all_results = {}
all_prot_seqs = {}
all_nucl_seqs = {}
for file in args.files:
    results = {}
    prot_seqs = {}
    nucl_seqs = {}
    # Is file nucleotide sequences?
    records = list(SeqIO.parse(file, 'fasta'))
    file_is_nucl = is_nucl(records[0])
    if file_is_nucl:
        nucl_records = records
        # Check if there is a corresponding protein file
        prot_file = f'{os.path.splitext(file)[0]}.faa'
        sys.stdout.write(f'{file} appears to be nucleotide sequences, checking for protein sequence file {prot_file}:\n')
        if os.path.exists(prot_file):
            sys.stdout.write(f'    Found.\nWill run hmmsearch on {prot_file}:\n')
            prot_records = list(SeqIO.parse(prot_file, 'fasta'))
        else:
            # If not, translate sequences and create stream of pretend file contents
            sys.stdout.write(f'    Not found. Will run hmmsearch on translated {file}:\n')
            prot_records = []
            for record in records: # Has to be a loop to avoid unknown ids
                prot_record = record.translate()
                prot_record.id = record.id
                prot_records.append(prot_record)
    else:
        # Check if there is a corresponding nucleotide file
        nucl_file = f'{os.path.splitext(file)[0]}.fna'
        sys.stdout.write(f'{file} appears to be protein sequences, checking for nucleotide sequence file {nucl_file}:\n')
        if os.path.exists(nucl_file):
            sys.stdout.write(f'    Found.\n')
            nucl_records = list(SeqIO.parse(nucl_file, 'fasta'))
        else:
            sys.stdout.write(f'    Not found. Will only output protein sequences.\n')
            nucl_records = None
        prot_records = list(SeqIO.parse(file, 'fasta'))

    # Create stream
    prot_stream = ''.join(f'>{x.id}\n{x.seq}\n' for x in prot_records)

    # Go through HMMs and save results
    for hmm in hmms.keys():
        sys.stdout.write(f'    {hmm}\n')
        results[hmm] = run_hmmsearch(hmms[hmm], prot_stream, cutoffs[hmm], args.threads, args.hmmsearch)
        if nucl_records is not None:
            if args.all:
                nucl_seqs[hmm] = [x for x in nucl_records if x.id in results[hmm].keys()]
            else:
                best_hit = max(results[hmm], key=results[hmm].get)
                nucl_seqs[hmm] = [x for x in nucl_records if x.id == best_hit]
        else:
            nucl_seqs[hmm] = None
        if args.all:
            prot_seqs[hmm] = [x for x in prot_records if x.id in results[hmm].keys()]
        else:
            best_hit = max(results[hmm], key=results[hmm].get)
            prot_seqs[hmm] = [x for x in prot_records if x.id == best_hit]

    all_results[file] = results
    all_prot_seqs[file] = prot_seqs
    all_nucl_seqs[file] = nucl_seqs

# Rearrange sequences
all_prot_seqs = {hmm:{file:all_prot_seqs[file][hmm] for file in all_prot_seqs.keys()} for hmm in hmms.keys()}
all_nucl_seqs = {hmm:{file:all_nucl_seqs[file][hmm] for file in all_nucl_seqs.keys()} for hmm in hmms.keys()}

# Output results table
os.makedirs(args.outdir, exist_ok=True)
with open(os.path.join(args.outdir, "marker_genes_scores.tsv"), 'w') as fo:
    for file in args.files:
        for hmm in hmms.keys():
            if args.all:
                for prot, score in all_results[file][hmm].items():
                    fo.write(f'{file}\t{hmm}\t{prot}\t{score}\n')
            else:
                best_hit = max(all_results[file][hmm], key=all_results[file][hmm].get)
                fo.write(f'{file}\t{hmm}\t{best_hit}\t{all_results[file][hmm][best_hit]}\n')

# Output sequences
for hmm in hmms.keys():
    with open(os.path.join(args.outdir, f'{hmm}.faa'), 'w') as fo:
        for file in args.files:
            filename = os.path.split(file)[1]
            for seq in all_prot_seqs[hmm][file]:
                fo.write(f'>{filename}_{seq.id}\n{seq.seq}\n')
