import pytest
import argparse
import tempfile
import io
import os
from fetchmgs import fetchMGs
from fetchmgs.test_config import *

def test_output_cutoffs_verybesthits():
    args = argparse.Namespace(
        v=True
    )
    hmms = {
        'COG0012': f'{PACKAGE_DIR}/data/COG0012.hmm',
        'COG0016': f'{PACKAGE_DIR}/data/COG0016.hmm',
        'COG0018': f'{PACKAGE_DIR}/data/COG0018.hmm',
        'COG0048': f'{PACKAGE_DIR}/data/COG0048.hmm',
        'COG0049': f'{PACKAGE_DIR}/data/COG0049.hmm',
        'COG0052': f'{PACKAGE_DIR}/data/COG0052.hmm',
        'COG0080': f'{PACKAGE_DIR}/data/COG0080.hmm',
        'COG0081': f'{PACKAGE_DIR}/data/COG0081.hmm',
        'COG0085': f'{PACKAGE_DIR}/data/COG0085.hmm',
        'COG0086': f'{PACKAGE_DIR}/data/COG0086.hmm',
        'COG0087': f'{PACKAGE_DIR}/data/COG0087.hmm',
        'COG0088': f'{PACKAGE_DIR}/data/COG0088.hmm',
        'COG0090': f'{PACKAGE_DIR}/data/COG0090.hmm',
        'COG0091': f'{PACKAGE_DIR}/data/COG0091.hmm',
        'COG0092': f'{PACKAGE_DIR}/data/COG0092.hmm',
        'COG0093': f'{PACKAGE_DIR}/data/COG0093.hmm',
        'COG0094': f'{PACKAGE_DIR}/data/COG0094.hmm',
        'COG0096': f'{PACKAGE_DIR}/data/COG0096.hmm',
        'COG0097': f'{PACKAGE_DIR}/data/COG0097.hmm',
        'COG0098': f'{PACKAGE_DIR}/data/COG0098.hmm',
        'COG0099': f'{PACKAGE_DIR}/data/COG0099.hmm',
        'COG0100': f'{PACKAGE_DIR}/data/COG0100.hmm',
        'COG0102': f'{PACKAGE_DIR}/data/COG0102.hmm',
        'COG0103': f'{PACKAGE_DIR}/data/COG0103.hmm',
        'COG0124': f'{PACKAGE_DIR}/data/COG0124.hmm',
        'COG0172': f'{PACKAGE_DIR}/data/COG0172.hmm',
        'COG0184': f'{PACKAGE_DIR}/data/COG0184.hmm',
        'COG0185': f'{PACKAGE_DIR}/data/COG0185.hmm',
        'COG0186': f'{PACKAGE_DIR}/data/COG0186.hmm',
        'COG0197': f'{PACKAGE_DIR}/data/COG0197.hmm',
        'COG0200': f'{PACKAGE_DIR}/data/COG0200.hmm',
        'COG0201': f'{PACKAGE_DIR}/data/COG0201.hmm',
        'COG0202': f'{PACKAGE_DIR}/data/COG0202.hmm',
        'COG0215': f'{PACKAGE_DIR}/data/COG0215.hmm',
        'COG0256': f'{PACKAGE_DIR}/data/COG0256.hmm',
        'COG0495': f'{PACKAGE_DIR}/data/COG0495.hmm',
        'COG0522': f'{PACKAGE_DIR}/data/COG0522.hmm',
        'COG0525': f'{PACKAGE_DIR}/data/COG0525.hmm',
        'COG0533': f'{PACKAGE_DIR}/data/COG0533.hmm',
        'COG0541': f'{PACKAGE_DIR}/data/COG0541.hmm',
        'COG0552': f'{PACKAGE_DIR}/data/COG0552.hmm'}
    new_cutoffs = {'COG0012': [510.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0016': [380.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0018': [530.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0048': [170.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0049': [210.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0052': [380.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0080': [230.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0081': [310.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0085': [1970.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0086': [870.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0087': [140.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0088': [240.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0090': [400.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0091': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0092': [360.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0093': [180.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0094': [220.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0096': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0097': [210.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0098': [220.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0099': [180.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0100': [180.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0102': [180.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0103': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0124': [460.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0172': [560.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0184': [110.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0185': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0186': [100.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0197': [190.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0200': [120.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0201': [440.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0202': [290.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0215': [580.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0256': [130.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0495': [1100.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0522': [250.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0525': [1270.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0533': [490.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0541': [570.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0552': [350.0, 5, 0, 0, 1.0, 1.0, 1.0]}

    with tempfile.TemporaryDirectory() as tmpdirname:
        args.o = tmpdirname
        fetchMGs.output_cutoffs(args=args, new_cutoffs=new_cutoffs, hmms=hmms)
        assert os.path.exists(os.path.join(tmpdirname, 'MG_BitScoreCutoffs.verybesthit.txt'))
        assert list(io.open(os.path.join(tmpdirname, 'MG_BitScoreCutoffs.verybesthit.txt'))) == list(io.open(f'{TRUE_OUTPUT_DIR}/cutoffs/MG_BitScoreCutoffs.verybesthit.txt'))

def test_output_cutoffs_allhits():
    args = argparse.Namespace(
        v=False
    )
    hmms = {
        'COG0012': f'{PACKAGE_DIR}/data/COG0012.hmm',
        'COG0016': f'{PACKAGE_DIR}/data/COG0016.hmm',
        'COG0018': f'{PACKAGE_DIR}/data/COG0018.hmm',
        'COG0048': f'{PACKAGE_DIR}/data/COG0048.hmm',
        'COG0049': f'{PACKAGE_DIR}/data/COG0049.hmm',
        'COG0052': f'{PACKAGE_DIR}/data/COG0052.hmm',
        'COG0080': f'{PACKAGE_DIR}/data/COG0080.hmm',
        'COG0081': f'{PACKAGE_DIR}/data/COG0081.hmm',
        'COG0085': f'{PACKAGE_DIR}/data/COG0085.hmm',
        'COG0086': f'{PACKAGE_DIR}/data/COG0086.hmm',
        'COG0087': f'{PACKAGE_DIR}/data/COG0087.hmm',
        'COG0088': f'{PACKAGE_DIR}/data/COG0088.hmm',
        'COG0090': f'{PACKAGE_DIR}/data/COG0090.hmm',
        'COG0091': f'{PACKAGE_DIR}/data/COG0091.hmm',
        'COG0092': f'{PACKAGE_DIR}/data/COG0092.hmm',
        'COG0093': f'{PACKAGE_DIR}/data/COG0093.hmm',
        'COG0094': f'{PACKAGE_DIR}/data/COG0094.hmm',
        'COG0096': f'{PACKAGE_DIR}/data/COG0096.hmm',
        'COG0097': f'{PACKAGE_DIR}/data/COG0097.hmm',
        'COG0098': f'{PACKAGE_DIR}/data/COG0098.hmm',
        'COG0099': f'{PACKAGE_DIR}/data/COG0099.hmm',
        'COG0100': f'{PACKAGE_DIR}/data/COG0100.hmm',
        'COG0102': f'{PACKAGE_DIR}/data/COG0102.hmm',
        'COG0103': f'{PACKAGE_DIR}/data/COG0103.hmm',
        'COG0124': f'{PACKAGE_DIR}/data/COG0124.hmm',
        'COG0172': f'{PACKAGE_DIR}/data/COG0172.hmm',
        'COG0184': f'{PACKAGE_DIR}/data/COG0184.hmm',
        'COG0185': f'{PACKAGE_DIR}/data/COG0185.hmm',
        'COG0186': f'{PACKAGE_DIR}/data/COG0186.hmm',
        'COG0197': f'{PACKAGE_DIR}/data/COG0197.hmm',
        'COG0200': f'{PACKAGE_DIR}/data/COG0200.hmm',
        'COG0201': f'{PACKAGE_DIR}/data/COG0201.hmm',
        'COG0202': f'{PACKAGE_DIR}/data/COG0202.hmm',
        'COG0215': f'{PACKAGE_DIR}/data/COG0215.hmm',
        'COG0256': f'{PACKAGE_DIR}/data/COG0256.hmm',
        'COG0495': f'{PACKAGE_DIR}/data/COG0495.hmm',
        'COG0522': f'{PACKAGE_DIR}/data/COG0522.hmm',
        'COG0525': f'{PACKAGE_DIR}/data/COG0525.hmm',
        'COG0533': f'{PACKAGE_DIR}/data/COG0533.hmm',
        'COG0541': f'{PACKAGE_DIR}/data/COG0541.hmm',
        'COG0552': f'{PACKAGE_DIR}/data/COG0552.hmm'}
    new_cutoffs = {'COG0012': [510.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0016': [380.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0018': [530.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0048': [170.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0049': [210.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0052': [380.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0080': [230.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0081': [310.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0085': [1970.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0086': [870.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0087': [140.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0088': [240.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0090': [400.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0091': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0092': [360.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0093': [180.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0094': [220.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0096': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0097': [210.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0098': [220.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0099': [180.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0100': [180.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0102': [180.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0103': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0124': [460.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0172': [560.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0184': [110.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0185': [150.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0186': [100.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0197': [190.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0200': [120.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0201': [440.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0202': [290.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0215': [580.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0256': [130.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0495': [1100.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0522': [250.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0525': [1270.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0533': [490.0, 5, 0, 0, 1.0, 1.0, 1.0], 'COG0541': [570.0, 5, 0, 0, 1.0, 1.0, 1.0],
     'COG0552': [350.0, 5, 0, 0, 1.0, 1.0, 1.0]}

    with tempfile.TemporaryDirectory() as tmpdirname:
        args.o = tmpdirname
        fetchMGs.output_cutoffs(args=args, new_cutoffs=new_cutoffs, hmms=hmms)
        assert os.path.exists(os.path.join(tmpdirname, 'MG_BitScoreCutoffs.allhits.txt'))
        assert list(io.open(os.path.join(tmpdirname, 'MG_BitScoreCutoffs.allhits.txt'))) == list(io.open(f'{TRUE_OUTPUT_DIR}/cutoffs/MG_BitScoreCutoffs.allhits.txt'))

