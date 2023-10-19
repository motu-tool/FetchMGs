from fetchmgs import fetchMGs
import unittest


class TestScoreCutoffs(unittest.TestCase):
    def test_score_cutoff(self):
        cutoff = 60.0
        neg = [61.0, 70.3, 73.1, 73.4, 77.0, 80.4]
        nvalid = 5
        pos = [516.8, 520.7, 541.3, 543.8, 562.1]

        assert fetchMGs.score_cutoff(pos, neg, nvalid, cutoff) == [60.0, 5, 6, 0, 0.45454545454545453, 1.0, 0.625]
            #     return ([cutoff, tp, fp, fn, precision, recall, fscore])

