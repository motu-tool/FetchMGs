import pytest

@pytest.fixture
def cutoffs():
    uncalibrated = {'COG0012': 60, 'COG0016': 60, 'COG0018': 60, 'COG0048': 60, 'COG0049': 60, 'COG0052': 60, 'COG0080': 60, 'COG0081': 60, 'COG0085': 60, 'COG0086': 60, 'COG0087': 60, 'COG0088': 60, 'COG0090': 60, 'COG0091': 60, 'COG0092': 60, 'COG0093': 60, 'COG0094': 60, 'COG0096': 60, 'COG0097': 60, 'COG0098': 60, 'COG0099': 60, 'COG0100': 60, 'COG0102': 60, 'COG0103': 60, 'COG0124': 60, 'COG0172': 60, 'COG0184': 60, 'COG0185': 60, 'COG0186': 60, 'COG0197': 60, 'COG0200': 60, 'COG0201': 60, 'COG0202': 60, 'COG0215': 60, 'COG0256': 60, 'COG0495': 60, 'COG0522': 60, 'COG0525': 60, 'COG0533': 60, 'COG0541': 60, 'COG0552': 60}
    allhits = {}
    besthits = {}
    return uncalibrated, allhits, besthits

@pytest.fixture
def hmms():
    return {'COG0012': 'fetchmgs/data/COG0012.hmm', 'COG0016': 'fetchmgs/data/COG0016.hmm', 'COG0018': 'fetchmgs/data/COG0018.hmm', 'COG0048': 'fetchmgs/data/COG0048.hmm', 'COG0049': 'fetchmgs/data/COG0049.hmm', 'COG0052': 'fetchmgs/data/COG0052.hmm', 'COG0080': 'fetchmgs/data/COG0080.hmm', 'COG0081': 'fetchmgs/data/COG0081.hmm', 'COG0085': 'fetchmgs/data/COG0085.hmm', 'COG0086': 'fetchmgs/data/COG0086.hmm', 'COG0087': 'fetchmgs/data/COG0087.hmm', 'COG0088': 'fetchmgs/data/COG0088.hmm', 'COG0090': 'fetchmgs/data/COG0090.hmm', 'COG0091': 'fetchmgs/data/COG0091.hmm', 'COG0092': 'fetchmgs/data/COG0092.hmm', 'COG0093': 'fetchmgs/data/COG0093.hmm', 'COG0094': 'fetchmgs/data/COG0094.hmm', 'COG0096': 'fetchmgs/data/COG0096.hmm', 'COG0097': 'fetchmgs/data/COG0097.hmm', 'COG0098': 'fetchmgs/data/COG0098.hmm', 'COG0099': 'fetchmgs/data/COG0099.hmm', 'COG0100': 'fetchmgs/data/COG0100.hmm', 'COG0102': 'fetchmgs/data/COG0102.hmm', 'COG0103': 'fetchmgs/data/COG0103.hmm', 'COG0124': 'fetchmgs/data/COG0124.hmm', 'COG0172': 'fetchmgs/data/COG0172.hmm', 'COG0184': 'fetchmgs/data/COG0184.hmm', 'COG0185': 'fetchmgs/data/COG0185.hmm', 'COG0186': 'fetchmgs/data/COG0186.hmm', 'COG0197': 'fetchmgs/data/COG0197.hmm', 'COG0200': 'fetchmgs/data/COG0200.hmm', 'COG0201': 'fetchmgs/data/COG0201.hmm', 'COG0202': 'fetchmgs/data/COG0202.hmm', 'COG0215': 'fetchmgs/data/COG0215.hmm', 'COG0256': 'fetchmgs/data/COG0256.hmm', 'COG0495': 'fetchmgs/data/COG0495.hmm', 'COG0522': 'fetchmgs/data/COG0522.hmm', 'COG0525': 'fetchmgs/data/COG0525.hmm', 'COG0533': 'fetchmgs/data/COG0533.hmm', 'COG0541': 'fetchmgs/data/COG0541.hmm', 'COG0552': 'fetchmgs/data/COG0552.hmm'}

