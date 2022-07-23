from . import Utils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''
FASTA file parser

---Attributes---
nt_id_nt_seq_records: Validated DNA sequence records with NCBI nucleotide accession numbers
nt_id_protein_seq_records: Validated protein sequence records with NCBI nucleotide accession numbers
protein_id_nt_seq_records: Validated DNA sequence records with NCBI protein accession numbers
protein_id_protein_seq_records: Validated protein sequence records with NCBI protein accession numbers
hybrid_seq_records: All validated DNA and protein records plus those without valid accession number,
                    for non-schema database file only
nt_seq_record_list: All validated DNA sequence records
protein_seq_record_list: All validated protein sequence records plus those without valid accession number
invalid_acc_num_fmt_seq_rec_ids: Sequence headers without valid accession number
unknown_seq_type_seq_rec_ids: Sequence headers for sequences having unknown sequence type
duplicated_headers: Duplicated sequence headers
redundant_seq_pairs: Sequence header pairs for identical sequences
seq_record_count: Total number of sequence records read
'''
class SequenceFileParser:
    def _init_data(self):
        self._nt_id_nt_seq_records = dict()
        self._nt_id_protein_seq_records = dict()
        self._protein_id_nt_seq_records = dict()
        self._protein_id_protein_seq_records = dict()
        self._hybrid_seq_records = dict()
        self._nt_seq_record_list = list()
        self._protein_seq_record_list = list()
        self._invalid_acc_num_fmt_seq_rec_ids = list()
        self._unknown_seq_type_seq_rec_ids = list()
        self._duplicated_headers = list()
        self._redundant_seq_pairs = list()
        self._seq_record_count = 0

    '''
    Function name: parse_seq_class_schema_db
    Inputs       : FASTA sequence file path, parse nucleotide sequence only flag
    Outputs      : Nil
    Description  : Parses schema ARG database file with validation for both sequence type
                   (DNA/protein/unknown) and accession number format (nucleotide/protein/missing)
    '''
    def parse_seq_class_schema_db(self, seq_file_path, is_nt_seq_only = False):
        self._init_data()

        with open(seq_file_path, 'r') as f:
            for seq_record in SeqIO.parse(f, 'fasta'):
                self._seq_record_count += 1
                seq_record.id = seq_record.description.replace(' ', '_')
                seq_record.seq = seq_record.seq.upper()

                '''Sequence type check'''
                seq_type = Utils.check_seq_type(str(seq_record.seq))
                if seq_type == '':
                    self._unknown_seq_type_seq_rec_ids.append(seq_record.description)
                    continue

                if is_nt_seq_only and seq_type != 'NT':
                    continue

                '''Accession number format check'''
                acc_num_type = Utils.check_acc_num_type(seq_record.id)
                if acc_num_type == '' and seq_type != 'AA':
                    self._invalid_acc_num_fmt_seq_rec_ids.append(seq_record.description)
                    continue

                if acc_num_type == 'NT':
                    nt_acc_num = Utils.extract_nt_acc_num(seq_record.id)
                    if seq_type == 'NT':
                        Utils.add_to_group(self._nt_id_nt_seq_records, nt_acc_num, seq_record)
                        self._nt_seq_record_list.append(seq_record)
                    elif seq_type == 'AA':
                        Utils.add_to_group(self._nt_id_protein_seq_records, nt_acc_num, seq_record)
                        self._protein_seq_record_list.append(seq_record)
                elif acc_num_type == 'AA':
                    protein_acc_num = Utils.extract_protein_acc_num(seq_record.id)
                    if seq_type == 'NT':
                        Utils.add_to_group(self._protein_id_nt_seq_records, protein_acc_num, seq_record)
                        self._nt_seq_record_list.append(seq_record)
                    elif seq_type == 'AA':
                        Utils.add_to_group(self._protein_id_protein_seq_records, protein_acc_num, seq_record)
                        self._protein_seq_record_list.append(seq_record)
                else:
                    '''
                    Sequence comparison/matching is performed at protein level, so only protein sequence
                    record list is updated
                    '''
                    if seq_type == 'AA':
                        self._protein_seq_record_list.append(seq_record)

    '''
    Function name: parse
    Inputs       : FASTA sequence file path, ignore accession number type flag
    Outputs      : Nil
    Description  : Parses non-schema ARG database file with validation for both sequence type
                   (DNA/protein/unknown) and accession number format (nucleotide/protein/missing)
    '''
    def parse(self, seq_file_path, is_ignore_acc_num_type = False):
        self._init_data()

        validated_seq_headers = set()
        validated_seq_strs = dict()

        with open(seq_file_path, 'r') as f:
            for seq_record in SeqIO.parse(f, 'fasta'):
                self._seq_record_count += 1
                seq_record.id = seq_record.description.replace(' ', '_')
                seq_record.seq = seq_record.seq.upper()

                '''Check for duplicated sequence header'''
                if Utils.is_duplicated_item(seq_record.id, validated_seq_headers):
                    self._duplicated_headers.append(seq_record.description)

                '''Sequence type check'''
                seq_type = Utils.check_seq_type(str(seq_record.seq))
                if seq_type == '':
                    self._unknown_seq_type_seq_rec_ids.append(seq_record.description)
                    continue
                else:
                    '''Check for redundant sequence, including reverse complementary sequence'''
                    redundant_seq_id, is_rev_comp = Utils.check_redundant_seq(seq_record.seq, seq_record.description,
                                                                              validated_seq_strs, seq_type)
                    if redundant_seq_id is not None:
                        self._redundant_seq_pairs.append((seq_record.description, redundant_seq_id, is_rev_comp))

                '''Accession number format check'''
                acc_num_type = Utils.check_acc_num_type(seq_record.id)
                if acc_num_type == '' and not is_ignore_acc_num_type:
                        self._invalid_acc_num_fmt_seq_rec_ids.append(seq_record.description)
                        continue

                if acc_num_type == 'NT':
                    nt_acc_num = Utils.extract_nt_acc_num(seq_record.id)
                    if seq_type == 'NT':
                        Utils.add_to_group(self._nt_id_nt_seq_records, nt_acc_num, seq_record)
                        Utils.add_to_group(self._hybrid_seq_records, nt_acc_num, seq_record)
                        self._nt_seq_record_list.append(seq_record)
                    elif seq_type == 'AA':
                        Utils.add_to_group(self._nt_id_protein_seq_records, nt_acc_num, seq_record)
                        Utils.add_to_group(self._hybrid_seq_records, nt_acc_num, seq_record)
                        self._protein_seq_record_list.append(seq_record)
                elif acc_num_type == 'AA':
                    protein_acc_num = Utils.extract_protein_acc_num(seq_record.id)
                    if seq_type == 'NT':
                        Utils.add_to_group(self._protein_id_nt_seq_records, protein_acc_num, seq_record)
                        Utils.add_to_group(self._hybrid_seq_records, protein_acc_num, seq_record)
                        self._nt_seq_record_list.append(seq_record)
                    elif seq_type == 'AA':
                        Utils.add_to_group(self._protein_id_protein_seq_records, protein_acc_num, seq_record)
                        Utils.add_to_group(self._hybrid_seq_records, protein_acc_num, seq_record)
                        self._protein_seq_record_list.append(seq_record)
                elif is_ignore_acc_num_type and acc_num_type == '':
                    '''
                    Sequence comparison/matching is performed at protein level, so only protein sequence
                    record list is updated
                    '''
                    if seq_type == 'NT':
                        Utils.add_to_group(self._hybrid_seq_records, seq_record.id, seq_record)
                    elif seq_type == 'AA':
                        Utils.add_to_group(self._hybrid_seq_records, seq_record.id, seq_record)
                        self._protein_seq_record_list.append(seq_record)

    def get_nt_id_nt_seq_records(self):
        return self._nt_id_nt_seq_records

    def get_nt_id_protein_seq_records(self):
        return self._nt_id_protein_seq_records

    def get_protein_id_nt_seq_records(self):
        return self._protein_id_nt_seq_records

    def get_protein_id_protein_seq_records(self):
        return self._protein_id_protein_seq_records

    def get_hybrid_seq_records(self):
        return self._hybrid_seq_records

    def get_nt_seq_record_list(self):
        return self._nt_seq_record_list

    def get_protein_seq_record_list(self):
        return self._protein_seq_record_list

    def get_invalid_acc_num_fmt_seq_rec_ids(self):
        return self._invalid_acc_num_fmt_seq_rec_ids

    def get_unknown_seq_type_seq_rec_ids(self):
        return self._unknown_seq_type_seq_rec_ids

    def get_duplicated_headers(self):
        return self._duplicated_headers

    def get_redundant_seq_pairs(self):
        return self._redundant_seq_pairs

    def get_seq_record_count(self):
        return self._seq_record_count
