import argparse
from Bio import SeqIO
import glob
import re
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='fetchMGs extracts the 40 single copy universal marker genes (decribed in Ciccarelli et al., Science, 2006 and Sorek et al., Science, 2007) from genomes and metagenomes in an easy and accurate manner.', formatter_class=argparse.RawTextHelpFormatter)

# Add -m option that does nothing for legacy
parser.add_argument('-m', '-mode', action='store_true', help='fetchMGs mode, see below')

subparsers = parser.add_subparsers(title='modes', description='valid modes', dest='mode')
subparsers.required = True

ext_parser = subparsers.add_parser('extraction', help='extract marker genes from sequences', formatter_class=argparse.RawTextHelpFormatter)
ext_parser.add_argument('file', help='multi-FASTA file with protein sequences from which universal single-copy marker genes should be extracted')
ext_parser.add_argument('-c', '-cog_used', nargs='+', default='all', help='orthologous group id(s) to be extracted; example: "COG0012"')
ext_parser.add_argument('-o', '-outdir', default='fetchmgs', help='output directory')
ext_parser.add_argument('-b', '-bitscore', default='data/MG_BitScoreCutoffs.allhits.txt', help='path to bitscore cutoff file')
ext_parser.add_argument('-l', '-library', default='data', help='path to directory that contains hmm models')
ext_parser.add_argument('-p', '-protein', action='store_true', help='set if nucleotide sequences file for <protein sequences> is not available')
ext_parser.add_argument('-d', '-dnaFastaFile', help='multi-FASTA file with nucleotide sequences; not neccesary if protein and nucleotide fasta file have the same name except .faa and .fna suffixes')
ext_parser.add_argument('-v', '-verybesthit_only', action='store_true', help='only extract the best hit of each COG from each genome\nrecommended to use, if extracting sequences from multiple reference genomes in the same file do not use it for metagenomes\nif this option is set fasta identifiers should be in the form: taxID.geneID and, if needed, have "project_id=XXX" in the header\nalternatively, set -i to ignore the headers; then, the best hit of each OG in the whole input file will be selected, regardless of the headers used')
ext_parser.add_argument('-i', '-ignore_headers', action='store_true', help='if this option is set in addition to -v, the best hit of each COG will be selected\nrecommended to use, if extracting sequences from a single genome in the same file')
ext_parser.add_argument('-t', '-threads', default=1, help='number of processors/threads to be used')
ext_parser.add_argument('-x', '-executable', default='bin', help='path to executables used by this script\ndefault is bin/\nif set to \'\' will search for executables in $PATH')

cal_parser = subparsers.add_parser('calibration', help='calibrate bitscores using results from extraction and a mapping file of known OGs')
cal_parser.add_argument('file', help='file with sequences that include marker genes (true positives)')
cal_parser.add_argument('map', help='tab-delimited file with true positive protein identifiers and COG ID')
cal_parser.add_argument('-o', '--outdir', default='fetchmgs', help='output directory')
cal_parser.add_argument('-b', '--bitscore', help='path to bitscore cutoff file')
# None of this mode implemented yet

# Arguments and modifications
args = parser.parse_args()
if args.x != '':
    if args.x[-1] != "/":
        args.x = f'{args.x}/'

# Compile list of HMM models
hmms = glob.glob(f'{args.l}/*.hmm')
if args.c == 'all':
    args.c = [os.path.splitext(os.path.split(x)[1])[0] for x in hmms]
args.c = sorted(args.c)
hmms = {x:f'{args.l}/{x}.hmm' for x in args.c}

# Parse cutoff file
if args.v:
    args.b = f'{args.l}/MG_BitScoreCutoffs.verybesthit.txt'
cutoffs = {}
with open(args.b, 'r') as fi:
    for line in fi.readlines():
        line = line.strip()
        if line[0] != '#':
           items = line.split('\t')
           cutoffs[items[0]] = int(items[1])

def run_hmmsearch(hmm_path, file, cutoff, threads, hmmsearch):
    '''
    Run hmmsearch from HMMER3, one hmm file against one set of protein records, then parse the results.
    '''
    try:
        output = subprocess.check_output(f'{hmmsearch}hmmsearch --noali --notextw --cpu {threads} -T {cutoff} {hmm_path} {file}', shell=True).decode()
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

# Organise files and records
prot_records = SeqIO.parse(args.file, 'fasta')
sys.stdout.write(f'Protein sequences: {args.file}\n')
# Check for nucleotide file
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

# Go through HMMs and save results
results = {}
for hmm in hmms.keys():
    sys.stdout.write(f'    {hmm}\n')
    results[hmm] = run_hmmsearch(hmms[hmm], args.file, cutoffs[hmm], args.t, args.x)
hit_ids = [k for result in results.values() for k in result.keys()]

# Get tax_ids if there are multiple genomes
if not args.i:
    if len(hit_ids[0].split(".")) > 1:
        tax_ids = set(x.split(".")[0] for x in hit_ids)
    else:
        tax_ids = None

# Filter if -v
if args.v:
    for hmm in hmms.keys():
        if tax_ids is not None:
            for tax_id in tax_ids:
                genome_results = {k:v for k,v in results[hmm].items() if k.split(".")[0]==tax_id}
                best_hit = max(genome_results, key=genome_results.get)
                [results[hmm].pop(k) for k in genome_results.keys() if k != best_hit]
        else:
            best_hit = max(results[hmm], key=results[hmm].get)
            results[hmm] = {best_hit:results[hmm][best_hit]}

# Only keep sequence records with hits
prot_hits = {x.id:x for x in prot_records if x.id in hit_ids}
if nucl_records is not None:
    nucl_hits = {x.id:x for x in nucl_records if x.id in hit_ids}

# Output results table
os.makedirs(args.o, exist_ok=True)
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
