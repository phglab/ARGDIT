from .Utils import check_seq_completeness

'''Object class for protein and its relevant information'''
class ProteinInfo:
    def __init__(self):
        self._acc_num = None
        self._name = None
        self._organism = None
        self._seq_str = None
        self._coding_gene_short_name = None
        self._is_5_partial = None
        self._is_3_partial = None

    @property
    def acc_num(self):
        return self._acc_num

    @acc_num.setter
    def acc_num(self, protein_acc_num):
        self._acc_num = protein_acc_num

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def organism(self):
        return self._organism

    @organism.setter
    def organism(self, organism):
        self._organism = organism

    @property
    def seq_str(self):
        return self._seq_str

    @seq_str.setter
    def seq_str(self, seq_str):
        if seq_str is not None:
            self._seq_str = seq_str.upper()

    @property
    def seq_len(self):
        return len(self._seq_str)

    @property
    def coding_gene_short_name(self):
        return self._coding_gene_short_name

    @coding_gene_short_name.setter
    def coding_gene_short_name(self, coding_gene_short_name):
        self._coding_gene_short_name = coding_gene_short_name

    @property
    def is_5_partial(self):
        return self._is_5_partial

    @is_5_partial.setter
    def is_5_partial(self, is_5_partial):
        self._is_5_partial = is_5_partial

    @property
    def is_3_partial(self):
        return self._is_3_partial

    @is_3_partial.setter
    def is_3_partial(self, is_3_partial):
        self._is_3_partial = is_3_partial

    def seq_completeness(self):
        return check_seq_completeness(self._is_5_partial, self._is_3_partial)
