from fetchmgs.fetchMGs import parse_cutoffs, import_files
from pathlib import Path
import pytest

TESTDIR = Path('/nfs/nas22/fs2202/biol_micro_bioinf/software/fetchMGs/example_datasets')

class FakeArgs:
    def __init__(self, params):
        self.mode = params.get('mode', None)
        self.file = params.get('file', None)
        self.c = params.get('c', None) # c
        self.l = params.get('l', None) # l
        self.b = params.get('b', None)
        self.v = params.get('v', None)
        self.d = params.get('d', None)
        self.map = params.get('map', None)
        self.p = params.get('p', None)

    def __repr__(self):
        return "\n".join([f"{k}: {v}" for k,v in zip(['mode', 'file', 'c', 'l', 'b', 'v', 'd', 'map'],
                                    [self.mode, self.file, self.c, self.l, self.b, self.v, self.d, self.map])])



fail_params = [(FakeArgs({
                          'mode': 'calibration',
                          'b': 'fetchmgs/data/MG_BitScoreCutoffs.allhits.txt',
                          }),
                "ERROR: When calibrating, the default bit score file provided with fetchMGs has to be used"),

                (FakeArgs({
                          'mode': 'extraction',
                          'b': 'fetchmgs/data/MG_BitScoreCutoffs.uncalibrated.txt',
                          'v': True
                          }),
                "ERROR: To extract using the option v|verybesthit, the file has to be calibrated with the -v option"),

                (FakeArgs({
                          'mode': 'extraction',
                          'b': 'fetchmgs/data/MG_BitScoreCutoffs.verybesthit.txt',
                          'v': False
                          }),
                "ERROR: To extract without the option v|verybesthit, the file has to be calibrated without the -v option")
              ]

@pytest.mark.parametrize("args, expected_error", fail_params)
def test_parse_cutoffs_fail(args, expected_error, capsys):
    with pytest.raises(SystemExit) as pytest_error:
        parse_cutoffs(args)
    out, err = capsys.readouterr()
    assert err == expected_error
    assert pytest_error.type == SystemExit


def test_parse_uncalibrated(cutoffs):
    args = FakeArgs(params= {'mode': 'calibration', 'b': 'fetchmgs/data/MG_BitScoreCutoffs.uncalibrated.txt'})
    output_cutoffs = parse_cutoffs(args)
    expected_cutoffs, _, _ = cutoffs
    assert output_cutoffs == expected_cutoffs

# todo same for  verybesthit and allhits (parametrize)

def test_import_files_defaults_calibration(cutoffs, hmms):
    args = FakeArgs({'l': 'fetchmgs/data',
                     'c': 'all',
                     'mode': 'calibration',
                     'file': TESTDIR/'example_data.faa',
                     'map':  'test/known_positives.map',
                      })
    out_hmms, out_cutoffs, valid_map, prot_records, nucl_records = import_files(args)
    expected_cutoffs, _, _ = cutoffs
    all_hmms = hmms
    assert out_cutoffs == expected_cutoffs
    assert out_hmms == all_hmms

# todo same with lots of different params