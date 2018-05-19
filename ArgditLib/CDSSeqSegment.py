'''Object class for CDS sequence segment'''
class CDSSeqSegment:
    def __init__(self, segment_seq_str, is_5_partial, is_3_partial, is_complementary):
        self._seq_str = segment_seq_str
        self._is_5_partial = is_5_partial
        self._is_3_partial = is_3_partial
        self._is_complementary = is_complementary

    @property
    def seq_str(self):
        return self._seq_str

    @property
    def length(self):
        return len(self._seq_str)

    @property
    def is_5_partial(self):
        return self._is_5_partial

    @property
    def is_3_partial(self):
        return self._is_3_partial

    @property
    def is_complementary(self):
        return self._is_complementary
