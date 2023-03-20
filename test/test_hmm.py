import os
import sys
import glob
import pytest
import tempfile
from pathlib import Path
from joblib import Parallel, delayed

# Directories required, modify if package structure changes
TEST_DIR       = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR     = f'{TEST_DIR}/../fetchmgs'         # Points to the package directory
EXAMPLE_DIR    = f'{TEST_DIR}/example_output/hmmResults/'      # HMM search example to be compared with tests
TEST_FILE      = f'{TEST_DIR}/test_data/example_data_genomes.faa'
TMP_DIR        = tempfile.TemporaryDirectory()

# Load fetchMGs functions
sys.path.insert(0, PARENT_DIR) 
from fetchMGs import run_hmmsearch, parse_cutoffs, parse_hmmsearch, extraction

# Locate the HMMSearch binary
HMMSEARCH_BIN  = f'{TEST_DIR}/../fetchmgs/bin/hmmsearch'    # Points to the hmmsearch bin in fetchmgs/bin
HMMSEARCH_BIN  = 'hmmsearch'    # Points to the hmmsearch bin in fetchmgs/bin

##########
# Dependencies
##########

class FakeArgs:
    """ Class required to interact with fetchMGs functions """
    def __init__(self, params):
        self.b    = params.get('b', None)            # Bitscore argument
        self.v    = params.get('v', None)            # Very best hit argument
        self.mode = params.get('mode', 'extraction') # Calibration or extraction, in this tests we only evaluate extraction 
        self.c    = params.get('c', 'all')           # COGs used, we consider the default 'all'

# Extraction modes requires cutoffs (2 different modes, depends on parse_cutoffs function) and hmms
cutoffs = {"allhits" : parse_cutoffs(FakeArgs({'b':f'{PARENT_DIR}/data/MG_BitScoreCutoffs.allhits.txt'})), 
           "besthits": parse_cutoffs(FakeArgs({'b':f'{PARENT_DIR}/data/MG_BitScoreCutoffs.verybesthit.txt', 'v':True}))}
hmms    = {cog_file.split('/')[-1].replace('.hmm', ''):cog_file for cog_file in glob.glob(f'{PARENT_DIR}/data/*.hmm')}   # This searches for all the hmm files in /data

# Functions to retrieve a fixture with test_results (2 different modes, depends on parse_cutoffs function) from hmms dictionary
def cmd_hmmsearh(cutoff, out_path, tbl_path, hmm_path):
    """ Command line call to run hmmsearch """
    cmd = f'{HMMSEARCH_BIN} --noali --notextw --cpu 1 -T {cutoff} -o {out_path} --domtblout {tbl_path} {hmm_path} {TEST_FILE}'
    os.system(cmd)
    return parse_hmmsearch(tbl_path)

def alternative_hmmsearch(v=0):
    """ 
    Run hmmsearch with the two set of cutoffs storing output files 
    in a tmpdir and return a dictionary of results. 
    Parallel implemented to speed up the testing process
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        arguments = []
        hmm_keys  = []
        for hmm, hmm_path in hmms.items():
            if v:
                cutoff = cutoffs['besthits'][hmm]
            else:
                cutoff = cutoffs['allhits'][hmm]
            arguments += [[cutoff, os.path.join(tmpdirname, f'{hmm}.out'), os.path.join(tmpdirname, f'{hmm}.dom'), hmm_path]]
            hmm_keys.append(hmm)
        results = Parallel(n_jobs=-1)(delayed(cmd_hmmsearh)(cutoff, out_path, tbl_path, hmm_path) for cutoff, out_path, tbl_path, hmm_path in arguments)
        return dict(zip(hmm_keys, results))

@pytest.fixture(scope = 'module', autouse=True)
def hmm_test_results():
    return  {'allhits':alternative_hmmsearch(0), 'besthits':alternative_hmmsearch(1)}

##########
# Actual tests
##########

def test_parse_hmmsearch():
    """
    Function to test the parsing of the standard output from hmmsearch.
    This function takes into account all *.dom files found in example_output.
    """
    hits = {}
    test_hits = {}
    for tbl_path in glob.glob(f'{TEST_DIR}/example_output/hmmResults/*.dom'):
        hits[tbl_path] = parse_hmmsearch(tbl_path)
        test_hits[tbl_path] = {}
        with open(tbl_path, 'r') as fi:
            for line in fi:
                if not line.startswith('#'):
                    items = line.strip().split()
                    test_hits[tbl_path][items[0]] = float(items[7])
    assert hits==test_hits, 'HMM search parsing test failed'

def test_run_hmmsearch_allhits(hmm_test_results):
    """ 
    Test hmmsearch works and it returns the same files for each hmm search when 
    using the allhits cutoffs.
    """  
    # Run hmmsearch with allhits
    results = {hmm:parse_hmmsearch(f'{TEST_DIR}/example_output/hmmResults/{hmm}.dom') for hmm in hmms.keys()}
    assert results==hmm_test_results['allhits'], f'HMM search using allhits cutoffs test failed, differences'

def test_run_hmmsearch_besthits(hmm_test_results):
    """ 
    Test hmmsearch works and it returns the same files for each hmm search when 
    using the besthits cutoffs.
    """  
    # Run hmmsearch with allhits
    results = {hmm:parse_hmmsearch(f'{TEST_DIR}/example_output_besthits/hmmResults/{hmm}.dom') for hmm in hmms.keys()}
    assert results==hmm_test_results['besthits'], f'HMM search using besthits cutoffs test failed'

def test_extraction_hit_ids_retrieval(hmm_test_results):
    """
    Test the retrieval of a list of hit ids (second output of the extraction function). 
    """
    results      = {hmm:parse_hmmsearch(f'{TEST_DIR}/example_output/hmmResults/{hmm}.dom') for hmm in hmms.keys()}
    hit_ids      = [k for result in results.values() for k in result.keys()]  # As done in the original function
    # We test only with one of the modes
    test_hit_ids = []
    for result in hmm_test_results['allhits'].values():
        for k in result.keys():
            test_hit_ids.append(k)
    assert hit_ids==test_hit_ids

def test_extraction_with_very_best_filter(hmm_test_results):
    """
    Test the retrieval hits when -v argument is passed. 
    This test is not ideal, to avoid running extraction function, 
    which would require to re-run all the hmmsearches, we copy the same
    function used in the original code to compare. 

    Ideally, the filtering should occur in a separate function.

    TODO: Multiple genomes testing
    """
    results      = {hmm:parse_hmmsearch(f'{TEST_DIR}/example_output/hmmResults/{hmm}.dom') for hmm in hmms.keys()}       # Precomputed
    test_results = {}
    for hmm in hmms.keys():
        # Original
        best_hit = max(results[hmm], key=results[hmm].get)
        results[hmm] = {best_hit:results[hmm][best_hit]}
        # Testing
        test_results[hmm] = {k:v for k, v in results[hmm].items() if v==max(results[hmm].values())}
    assert results==test_results