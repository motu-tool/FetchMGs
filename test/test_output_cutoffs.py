import unittest
import pathlib
import argparse
import io
import os
from fetchmgs import fetchMGs

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PKG_DIR  = f'{TEST_DIR}/../fetchmgs'

class TestOutputCutoffs(unittest.TestCase):
    def test_output_cutoffs_verybesthits(self):
        mock_args = argparse.Namespace(
            o=pathlib.Path('.'),
            v=True
        )
        hmms = {
            'COG0012': pathlib.Path(f'{PKG_DIR}/data/COG0012.hmm'),
            'COG0016': pathlib.Path(f'{PKG_DIR}/data/COG0016.hmm'),
            'COG0018': pathlib.Path(f'{PKG_DIR}/data/COG0018.hmm'),
            'COG0048': pathlib.Path(f'{PKG_DIR}/data/COG0048.hmm'),
            'COG0049': pathlib.Path(f'{PKG_DIR}/data/COG0049.hmm'),
            'COG0052': pathlib.Path(f'{PKG_DIR}/data/COG0052.hmm'),
            'COG0080': pathlib.Path(f'{PKG_DIR}/data/COG0080.hmm'),
            'COG0081': pathlib.Path(f'{PKG_DIR}/data/COG0081.hmm'),
            'COG0085': pathlib.Path(f'{PKG_DIR}/data/COG0085.hmm'),
            'COG0086': pathlib.Path(f'{PKG_DIR}/data/COG0086.hmm'),
            'COG0087': pathlib.Path(f'{PKG_DIR}/data/COG0087.hmm'),
            'COG0088': pathlib.Path(f'{PKG_DIR}/data/COG0088.hmm'),
            'COG0090': pathlib.Path(f'{PKG_DIR}/data/COG0090.hmm'),
            'COG0091': pathlib.Path(f'{PKG_DIR}/data/COG0091.hmm'),
            'COG0092': pathlib.Path(f'{PKG_DIR}/data/COG0092.hmm'),
            'COG0093': pathlib.Path(f'{PKG_DIR}/data/COG0093.hmm'),
            'COG0094': pathlib.Path(f'{PKG_DIR}/data/COG0094.hmm'),
            'COG0096': pathlib.Path(f'{PKG_DIR}/data/COG0096.hmm'),
            'COG0097': pathlib.Path(f'{PKG_DIR}/data/COG0097.hmm'),
            'COG0098': pathlib.Path(f'{PKG_DIR}/data/COG0098.hmm'),
            'COG0099': pathlib.Path(f'{PKG_DIR}/data/COG0099.hmm'),
            'COG0100': pathlib.Path(f'{PKG_DIR}/data/COG0100.hmm'),
            'COG0102': pathlib.Path(f'{PKG_DIR}/data/COG0102.hmm'),
            'COG0103': pathlib.Path(f'{PKG_DIR}/data/COG0103.hmm'),
            'COG0124': pathlib.Path(f'{PKG_DIR}/data/COG0124.hmm'),
            'COG0172': pathlib.Path(f'{PKG_DIR}/data/COG0172.hmm'),
            'COG0184': pathlib.Path(f'{PKG_DIR}/data/COG0184.hmm'),
            'COG0185': pathlib.Path(f'{PKG_DIR}/data/COG0185.hmm'),
            'COG0186': pathlib.Path(f'{PKG_DIR}/data/COG0186.hmm'),
            'COG0197': pathlib.Path(f'{PKG_DIR}/data/COG0197.hmm'),
            'COG0200': pathlib.Path(f'{PKG_DIR}/data/COG0200.hmm'),
            'COG0201': pathlib.Path(f'{PKG_DIR}/data/COG0201.hmm'),
            'COG0202': pathlib.Path(f'{PKG_DIR}/data/COG0202.hmm'),
            'COG0215': pathlib.Path(f'{PKG_DIR}/data/COG0215.hmm'),
            'COG0256': pathlib.Path(f'{PKG_DIR}/data/COG0256.hmm'),
            'COG0495': pathlib.Path(f'{PKG_DIR}/data/COG0495.hmm'),
            'COG0522': pathlib.Path(f'{PKG_DIR}/data/COG0522.hmm'),
            'COG0525': pathlib.Path(f'{PKG_DIR}/data/COG0525.hmm'),
            'COG0533': pathlib.Path(f'{PKG_DIR}/data/COG0533.hmm'),
            'COG0541': pathlib.Path(f'{PKG_DIR}/data/COG0541.hmm'),
            'COG0552': pathlib.Path(f'{PKG_DIR}/data/COG0552.hmm')}
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

        fetchMGs.output_cutoffs(args=mock_args, new_cutoffs=new_cutoffs, hmms=hmms)
        self.assertTrue(pathlib.Path('MG_BitScoreCutoffs.verybesthit.txt').exists())
        self.assertListEqual(
            list(io.open('MG_BitScoreCutoffs.verybesthit.txt')),
            list(io.open(pathlib.Path(f'{TEST_DIR}/output/MG_BitScoreCutoffs.verybesthit.txt'))))
        pathlib.Path('MG_BitScoreCutoffs.verybesthit.txt').unlink()

    def test_output_cutoffs_allhits(self):
        mock_args = argparse.Namespace(
            o=pathlib.Path('.'),
            v=False
        )
        hmms = {
            'COG0012': pathlib.Path(f'{PKG_DIR}/data/COG0012.hmm'),
            'COG0016': pathlib.Path(f'{PKG_DIR}/data/COG0016.hmm'),
            'COG0018': pathlib.Path(f'{PKG_DIR}/data/COG0018.hmm'),
            'COG0048': pathlib.Path(f'{PKG_DIR}/data/COG0048.hmm'),
            'COG0049': pathlib.Path(f'{PKG_DIR}/data/COG0049.hmm'),
            'COG0052': pathlib.Path(f'{PKG_DIR}/data/COG0052.hmm'),
            'COG0080': pathlib.Path(f'{PKG_DIR}/data/COG0080.hmm'),
            'COG0081': pathlib.Path(f'{PKG_DIR}/data/COG0081.hmm'),
            'COG0085': pathlib.Path(f'{PKG_DIR}/data/COG0085.hmm'),
            'COG0086': pathlib.Path(f'{PKG_DIR}/data/COG0086.hmm'),
            'COG0087': pathlib.Path(f'{PKG_DIR}/data/COG0087.hmm'),
            'COG0088': pathlib.Path(f'{PKG_DIR}/data/COG0088.hmm'),
            'COG0090': pathlib.Path(f'{PKG_DIR}/data/COG0090.hmm'),
            'COG0091': pathlib.Path(f'{PKG_DIR}/data/COG0091.hmm'),
            'COG0092': pathlib.Path(f'{PKG_DIR}/data/COG0092.hmm'),
            'COG0093': pathlib.Path(f'{PKG_DIR}/data/COG0093.hmm'),
            'COG0094': pathlib.Path(f'{PKG_DIR}/data/COG0094.hmm'),
            'COG0096': pathlib.Path(f'{PKG_DIR}/data/COG0096.hmm'),
            'COG0097': pathlib.Path(f'{PKG_DIR}/data/COG0097.hmm'),
            'COG0098': pathlib.Path(f'{PKG_DIR}/data/COG0098.hmm'),
            'COG0099': pathlib.Path(f'{PKG_DIR}/data/COG0099.hmm'),
            'COG0100': pathlib.Path(f'{PKG_DIR}/data/COG0100.hmm'),
            'COG0102': pathlib.Path(f'{PKG_DIR}/data/COG0102.hmm'),
            'COG0103': pathlib.Path(f'{PKG_DIR}/data/COG0103.hmm'),
            'COG0124': pathlib.Path(f'{PKG_DIR}/data/COG0124.hmm'),
            'COG0172': pathlib.Path(f'{PKG_DIR}/data/COG0172.hmm'),
            'COG0184': pathlib.Path(f'{PKG_DIR}/data/COG0184.hmm'),
            'COG0185': pathlib.Path(f'{PKG_DIR}/data/COG0185.hmm'),
            'COG0186': pathlib.Path(f'{PKG_DIR}/data/COG0186.hmm'),
            'COG0197': pathlib.Path(f'{PKG_DIR}/data/COG0197.hmm'),
            'COG0200': pathlib.Path(f'{PKG_DIR}/data/COG0200.hmm'),
            'COG0201': pathlib.Path(f'{PKG_DIR}/data/COG0201.hmm'),
            'COG0202': pathlib.Path(f'{PKG_DIR}/data/COG0202.hmm'),
            'COG0215': pathlib.Path(f'{PKG_DIR}/data/COG0215.hmm'),
            'COG0256': pathlib.Path(f'{PKG_DIR}/data/COG0256.hmm'),
            'COG0495': pathlib.Path(f'{PKG_DIR}/data/COG0495.hmm'),
            'COG0522': pathlib.Path(f'{PKG_DIR}/data/COG0522.hmm'),
            'COG0525': pathlib.Path(f'{PKG_DIR}/data/COG0525.hmm'),
            'COG0533': pathlib.Path(f'{PKG_DIR}/data/COG0533.hmm'),
            'COG0541': pathlib.Path(f'{PKG_DIR}/data/COG0541.hmm'),
            'COG0552': pathlib.Path(f'{PKG_DIR}/data/COG0552.hmm')}
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

        fetchMGs.output_cutoffs(args=mock_args, new_cutoffs=new_cutoffs, hmms=hmms)
        self.assertTrue(pathlib.Path('MG_BitScoreCutoffs.allhits.txt').exists())
        self.assertListEqual(
            list(io.open('MG_BitScoreCutoffs.allhits.txt')),
            list(io.open(pathlib.Path(f'{TEST_DIR}/output/MG_BitScoreCutoffs.allhits.txt'))))
        pathlib.Path('MG_BitScoreCutoffs.allhits.txt').unlink()


if __name__ == '__main__':
    unittest.main()
