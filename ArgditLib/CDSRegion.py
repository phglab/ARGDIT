from .Utils import check_seq_completeness

'''Object class for CDS region downloaded from NCBI nucleotide database'''
class CDSRegion:
    def __init__(self):
        self._region_ranges = list()
        self._length = 0
        self._codon_start = 1
        self._protein_id = None
        self._gene_short_name = None
        self._product = None
        self._is_5_partial = False
        self._is_3_partial = False
        self._is_enforce_start_codon = False

    def add_region_range(self, range_start, range_end):
        self._region_ranges.append((range_start, range_end))
        self._length += (abs(range_end - range_start) + 1)

    @property
    def codon_start(self):
        return self._codon_start

    @codon_start.setter
    def codon_start(self, codon_start):
        self._codon_start = codon_start

    @property
    def protein_id(self):
        return self._protein_id

    @protein_id.setter
    def protein_id(self, protein_acc_num):
        self._protein_id = protein_acc_num

    @property
    def gene_short_name(self):
        return self._gene_short_name

    @gene_short_name.setter
    def gene_short_name(self, gene_short_name):
        self._gene_short_name = gene_short_name

    @property
    def product(self):
        return self._product

    @product.setter
    def product(self, product):
        self._product = product

    @property
    def is_5_partial(self):
        return self._is_5_partial

    @is_5_partial.setter
    def is_5_partial(self, is_5_partial_cds):
        self._is_5_partial = is_5_partial_cds

    @property
    def is_3_partial(self):
        return self._is_3_partial

    @is_3_partial.setter
    def is_3_partial(self, is_3_partial_cds):
        self._is_3_partial = is_3_partial_cds

    def seq_completeness(self):
        return check_seq_completeness(self._is_5_partial, self.is_3_partial)

    @property
    def is_enforce_start_codon(self):
        return self._is_enforce_start_codon

    def set_enforce_start_codon(self, translate_except_codon_start):
        ordered_region_ranges = sorted(self._region_ranges, key = lambda region_range:region_range[0])
        if self.is_complementary:
            self._is_enforce_start_codon = (ordered_region_ranges[-1][1] == translate_except_codon_start)
        else:
            self._is_enforce_start_codon = (ordered_region_ranges[0][0] == translate_except_codon_start)

    @property
    def length(self):
        return self._length

    @property
    def codon_count(self):
        return (self._length - (self._codon_start - 1)) // 3

    @property
    def is_complementary(self):
        if self._length > 0:
            return self._region_ranges[0][1] < self._region_ranges[0][0]
        else:
            return False

    def _sort_region_ranges(self):
        if self.is_complementary:
            return sorted(self._region_ranges, key = lambda region_range:region_range[0], reverse = True)
        else:
            return sorted(self._region_ranges, key = lambda region_range:region_range[0])

    def get_region_range(self):
        ordered_region_ranges = self._sort_region_ranges()

        return ordered_region_ranges[0][0], ordered_region_ranges[-1][1]

    def get_region_range_annotation(self):
        ordered_region_ranges = self._sort_region_ranges()
        region_range_annotation = None

        for region_range in ordered_region_ranges:
            if region_range_annotation is None:
                region_range_annotation = '{}-{}'.format(region_range[0], region_range[1])
            else:
                region_range_annotation = '{},{}-{}'.format(region_range_annotation, region_range[0], region_range[1])
        if region_range_annotation is None:
            return ''
        else:
            return region_range_annotation
