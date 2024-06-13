from fetchmgs import fetchmgs
from fetchmgs.test_config import *
import os
import pytest

class FakeArgs:
    """
    Class to mimic args.parse object (i.e. parsed params passed by the user)

    """

    def __init__(self, params):
        self.mode = params.get("mode", None)
        self.file = params.get("file", None)
        self.c = params.get("c", None)
        self.l = params.get("l", None)
        self.b = params.get("b", None)
        self.v = params.get("v", None)
        self.d = params.get("d", None)
        self.map = params.get("map", None)
        self.p = params.get("p", None)
        self.o = params.get("o", "output")

    def __repr__(self):
        return "\n".join(
            [
                f"{k}: {v}"
                for k, v in zip(
                    ["mode", "file", "c", "l", "b", "v", "d", "map", "o"],
                    [
                        self.mode,
                        self.file,
                        self.c,
                        self.l,
                        self.b,
                        self.v,
                        self.d,
                        self.map,
                        self.o
                    ],
                )
            ]
        )


# Params that  should lead to error

fail_params = [
    (
        FakeArgs(
            {
                "mode": "calibration",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.allhits.txt',
            }
        ),
        "ERROR: When calibrating, the default bit score file provided with FetchMGs has to be used",
    ),
    (
        FakeArgs(
            {
                "mode": "extraction",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.uncalibrated.txt',
                "v": True,
            }
        ),
        "ERROR: To extract using the option v|verybesthit, the file has to be calibrated with the -v option",
    ),
    (
        FakeArgs(
            {
                "mode": "extraction",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.verybesthit.txt',
                "v": False,
            }
        ),
        "ERROR: To extract without the option v|verybesthit, the file has to be calibrated without the -v option",
    ),
]


@pytest.mark.parametrize("args, expected_error", fail_params)
def test_parse_cutoffs_fail(args, expected_error, capsys):
    with pytest.raises(SystemExit) as pytest_error:
        fetchmgs.parse_cutoffs(args)
    out, err = capsys.readouterr()
    assert err == expected_error
    assert pytest_error.type == SystemExit


# Define expected cutoffs and locations of  files.
cutoffs = {
    "uncalibrated": {
        "COG0012": 60,
        "COG0016": 60,
        "COG0018": 60,
        "COG0048": 60,
        "COG0049": 60,
        "COG0052": 60,
        "COG0080": 60,
        "COG0081": 60,
        "COG0085": 60,
        "COG0086": 60,
        "COG0087": 60,
        "COG0088": 60,
        "COG0090": 60,
        "COG0091": 60,
        "COG0092": 60,
        "COG0093": 60,
        "COG0094": 60,
        "COG0096": 60,
        "COG0097": 60,
        "COG0098": 60,
        "COG0099": 60,
        "COG0100": 60,
        "COG0102": 60,
        "COG0103": 60,
        "COG0124": 60,
        "COG0172": 60,
        "COG0184": 60,
        "COG0185": 60,
        "COG0186": 60,
        "COG0197": 60,
        "COG0200": 60,
        "COG0201": 60,
        "COG0202": 60,
        "COG0215": 60,
        "COG0256": 60,
        "COG0495": 60,
        "COG0522": 60,
        "COG0525": 60,
        "COG0533": 60,
        "COG0541": 60,
        "COG0552": 60,
    },
    "allhits": {
        "COG0012": 210,
        "COG0016": 240,
        "COG0018": 340,
        "COG0048": 100,
        "COG0049": 120,
        "COG0052": 140,
        "COG0080": 90,
        "COG0081": 130,
        "COG0085": 1020,
        "COG0086": 60,
        "COG0087": 120,
        "COG0088": 110,
        "COG0090": 180,
        "COG0091": 80,
        "COG0092": 120,
        "COG0093": 80,
        "COG0094": 110,
        "COG0096": 80,
        "COG0097": 100,
        "COG0098": 140,
        "COG0099": 120,
        "COG0100": 80,
        "COG0102": 100,
        "COG0103": 80,
        "COG0124": 320,
        "COG0172": 170,
        "COG0184": 60,
        "COG0185": 70,
        "COG0186": 80,
        "COG0197": 70,
        "COG0200": 60,
        "COG0201": 210,
        "COG0202": 80,
        "COG0215": 400,
        "COG0256": 70,
        "COG0495": 450,
        "COG0522": 80,
        "COG0525": 740,
        "COG0533": 300,
        "COG0541": 450,
        "COG0552": 220,
    },
    "besthits": {
        "COG0012": 210,
        "COG0016": 240,
        "COG0018": 340,
        "COG0048": 100,
        "COG0049": 120,
        "COG0052": 140,
        "COG0080": 90,
        "COG0081": 130,
        "COG0085": 540,
        "COG0086": 60,
        "COG0087": 120,
        "COG0088": 110,
        "COG0090": 180,
        "COG0091": 60,
        "COG0092": 120,
        "COG0093": 80,
        "COG0094": 110,
        "COG0096": 80,
        "COG0097": 100,
        "COG0098": 140,
        "COG0099": 120,
        "COG0100": 80,
        "COG0102": 100,
        "COG0103": 80,
        "COG0124": 220,
        "COG0172": 170,
        "COG0184": 60,
        "COG0185": 70,
        "COG0186": 80,
        "COG0197": 70,
        "COG0200": 60,
        "COG0201": 210,
        "COG0202": 80,
        "COG0215": 380,
        "COG0256": 70,
        "COG0495": 420,
        "COG0522": 80,
        "COG0525": 740,
        "COG0533": 300,
        "COG0541": 450,
        "COG0552": 220,
    },
}
hmms = {
    "COG0012": "fetchmgs/data/COG0012.hmm",
    "COG0016": "fetchmgs/data/COG0016.hmm",
    "COG0018": "fetchmgs/data/COG0018.hmm",
    "COG0048": "fetchmgs/data/COG0048.hmm",
    "COG0049": "fetchmgs/data/COG0049.hmm",
    "COG0052": "fetchmgs/data/COG0052.hmm",
    "COG0080": "fetchmgs/data/COG0080.hmm",
    "COG0081": "fetchmgs/data/COG0081.hmm",
    "COG0085": "fetchmgs/data/COG0085.hmm",
    "COG0086": "fetchmgs/data/COG0086.hmm",
    "COG0087": "fetchmgs/data/COG0087.hmm",
    "COG0088": "fetchmgs/data/COG0088.hmm",
    "COG0090": "fetchmgs/data/COG0090.hmm",
    "COG0091": "fetchmgs/data/COG0091.hmm",
    "COG0092": "fetchmgs/data/COG0092.hmm",
    "COG0093": "fetchmgs/data/COG0093.hmm",
    "COG0094": "fetchmgs/data/COG0094.hmm",
    "COG0096": "fetchmgs/data/COG0096.hmm",
    "COG0097": "fetchmgs/data/COG0097.hmm",
    "COG0098": "fetchmgs/data/COG0098.hmm",
    "COG0099": "fetchmgs/data/COG0099.hmm",
    "COG0100": "fetchmgs/data/COG0100.hmm",
    "COG0102": "fetchmgs/data/COG0102.hmm",
    "COG0103": "fetchmgs/data/COG0103.hmm",
    "COG0124": "fetchmgs/data/COG0124.hmm",
    "COG0172": "fetchmgs/data/COG0172.hmm",
    "COG0184": "fetchmgs/data/COG0184.hmm",
    "COG0185": "fetchmgs/data/COG0185.hmm",
    "COG0186": "fetchmgs/data/COG0186.hmm",
    "COG0197": "fetchmgs/data/COG0197.hmm",
    "COG0200": "fetchmgs/data/COG0200.hmm",
    "COG0201": "fetchmgs/data/COG0201.hmm",
    "COG0202": "fetchmgs/data/COG0202.hmm",
    "COG0215": "fetchmgs/data/COG0215.hmm",
    "COG0256": "fetchmgs/data/COG0256.hmm",
    "COG0495": "fetchmgs/data/COG0495.hmm",
    "COG0522": "fetchmgs/data/COG0522.hmm",
    "COG0525": "fetchmgs/data/COG0525.hmm",
    "COG0533": "fetchmgs/data/COG0533.hmm",
    "COG0541": "fetchmgs/data/COG0541.hmm",
    "COG0552": "fetchmgs/data/COG0552.hmm",
}
first_sequence = {
    "seq_id": "gene_1",
    "nt_seq": "ATGTTACCAAGGGACGATTCTGAAAGCCCATTAGCTTTACAATTGTTTGGTGGCAATAAGGACATAATATTAAATGCTATAGATTTTGTTGAAAAAGAAGCTAAATATGATTTCTTAGATTTTAATATGGGTTGTCCTGTTCCTAAAGTTATGAAACAAGAGGCTGGTTCGTATTGGTTAAAACGACAAGATGAAGTATATGATTTGCTTCATTCAATGGTACAAAAATCTCATAAGCCTGTAATAATTAAAATTAGACTTGGATTTGATAAAAAGCATATTAATGCAGTAGAAATTGCTAAAATTGCTGAACAAGCTGGAATTAAAGCTTTAGCTGTTCATGGAAGAACAAGAGATGAATACTATAACGGTTTTCCACATTATGAGGAAATTGCTAAAGTTAAAAATAGTGTTTCAATACCGGTTATAGCCAATGGAAATATTGATTTAAACAATATTGCTGAAGTTGAAAAAATAACGAATGCTGATGCATTTATGATTGGCAGAGGATGTCTTGGTAACCCACTAATATTTACTGATCTTATTAATTCGGAAGAAGGTAAACCGTTTATAGAACACACATTAGAAAAACAAATAAATTTGATGAAAAAACATTTTGATATGATTATAAATTATATGGGTGAATATAACGGGGTAAGATACTTTAGAGGCTTATCAGTTTTATATGTTAAGGGATTTGACAACGCAAAATATATTAAAAATAAATTGATATCAATGAATACAGTTGATGATTTTAACTCAATAATAAAAGAGACTAAAGATTACTATTTGAAATAA",
    "aa_seq": "MLPRDDSESPLALQLFGGNKDIILNAIDFVEKEAKYDFLDFNMGCPVPKVMKQEAGSYWLKRQDEVYDLLHSMVQKSHKPVIIKIRLGFDKKHINAVEIAKIAEQAGIKALAVHGRTRDEYYNGFPHYEEIAKVKNSVSIPVIANGNIDLNNIAEVEKITNADAFMIGRGCLGNPLIFTDLINSEEGKPFIEHTLEKQINLMKKHFDMIINYMGEYNGVRYFRGLSVLYVKGFDNAKYIKNKLISMNTVDDFNSIIKETKDYYLK",
}


def pass_cutoff_params():
    """
    Generates parameters and expected outcomes for cuoff_params function
    """
    uncal_args, allhits_args, besthits_args = [
        FakeArgs(
            params={
                "mode": "calibration",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.uncalibrated.txt',
            }
        ),
        FakeArgs(
            params={
                "mode": "extraction",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.allhits.txt',
                "v": False,
            }
        ),
        FakeArgs(
            params={
                "mode": "extraction",
                "b": f'{PACKAGE_DIR}/data/MG_BitScoreCutoffs.verybesthit.txt',
                "v": True,
            }
        ),
    ]
    return [
        (uncal_args, cutoffs["uncalibrated"]),
        (allhits_args, cutoffs["allhits"]),
        (besthits_args, cutoffs["besthits"]),
    ]


# Testing parse_cutoffs
@pytest.mark.parametrize("args, expected_cutoffs", pass_cutoff_params())
def test_parse_cutoffs(args, expected_cutoffs):
    output_cutoffs = fetchmgs.parse_cutoffs(args)
    assert output_cutoffs == expected_cutoffs


def pass_setup_params():
    """
    Generates params and expected outcomes for setup function
    """
    cal_args, besthit_args, allhit_args = [
        FakeArgs(
            {
                "l": "fetchmgs/data",
                "c": "all",
                "mode": "calibration",
                "file": f'{TEST_INPUT_DIR}/setup_test_data.faa',
            }
        ),
        FakeArgs(
            {
                "l": "fetchmgs/data",
                "c": "all",
                "mode": "extraction",
                "v": True,
                "file": f'{TEST_INPUT_DIR}/setup_test_data.faa',
            }
        ),
        FakeArgs(
            {
                "l": "fetchmgs/data",
                "c": "all",
                "mode": "extraction",
                "v": False,
                "file": f'{TEST_INPUT_DIR}/setup_test_data.faa',
            }
        ),
    ]
    return [
        (cal_args, cutoffs["uncalibrated"], hmms),
        (besthit_args, cutoffs["besthits"], hmms),
        (allhit_args, cutoffs["allhits"], hmms),
    ]


# Test setup
@pytest.mark.parametrize(
    "args, expected_cutoffs, expected_hmms", pass_setup_params()
)
def test_setup(args, expected_cutoffs, expected_hmms):
    out_hmms, out_cutoffs, prot_records, nucl_records = fetchmgs.setup(args)
    fseq = next(prot_records)
    assert out_cutoffs == expected_cutoffs
    assert out_hmms == expected_hmms
    assert fseq.id == first_sequence["seq_id"]
    assert fseq.seq == first_sequence["aa_seq"]


# Todos:
# todo test c != 'all'
# todo check seq in all setup tests
