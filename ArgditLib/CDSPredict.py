'''Module to predict CDS sequence segment from input nucleotide/protein sequence'''

from Bio.Seq import Seq
from .Translate import Translate
from .CDSSeqSegment import CDSSeqSegment
from .Utils import trim_version, add_to_group

'''
Function name: _detect_start_and_stop_codons
Inputs       : Nucleotide sequence string
Outputs      : All detected start and stop codons
Description  : Detects all start codons within the first 5 nucleotides, and all stop codons within
               the last 5 nucleotides
'''
def _detect_start_and_stop_codons(seq_str):
    cds_starts = list()

    for i in range(3):
        if seq_str[i:i + 3] in Translate.get_start_codons():
            cds_starts.append(i)

    cds_ends = list()
    seq_len = len(seq_str)

    for i in range(seq_len - 5, seq_len - 2):
        if seq_str[i:i + 3] in Translate.get_stop_codons():
            cds_ends.append(i + 3)

    return cds_starts, cds_ends

'''
Function name: _generate_cds_seq_segments
Inputs       : Nucleotide sequence string, sequence complementarity
Outputs      : Candidate CDS sequence segments
Description  : Generates candidate CDS sequence segments according to four cases of CDS completeness:
               1. Complete
               2. 5' partial
               3. 3' partial
               4. 5' and 3' partial
'''
def _generate_cds_seq_segments(seq_str, is_complementary):
    cds_starts, cds_ends = _detect_start_and_stop_codons(seq_str)

    cds_seq_segments = list()
    '''Assume the CDS is 5' and 3' partial'''
    cds_seq_segments.append(CDSSeqSegment(seq_str, True, True, is_complementary))

    seq_len = len(seq_str)

    for cds_start in cds_starts:
        for cds_end in cds_ends:
            '''When the first and last codons correspond to start and stop codon, respectively, and are
               separated by N codons, create a complete CDS'''
            if (cds_end - cds_start) % 3 == 0:
                cds_seq_segments.append(CDSSeqSegment(seq_str[cds_start:cds_end], False, False, is_complementary))

        '''For each start codon, create a 3' partial CDS'''
        cds_seq_segments.append(CDSSeqSegment(seq_str[cds_start:seq_len], False, True, is_complementary))

    for cds_end in cds_ends:
        '''For each stop codon, create a 5' partial CDS'''
        cds_seq_segments.append(CDSSeqSegment(seq_str[0:cds_end], True, False, is_complementary))

    return cds_seq_segments

'''
Function name: _predict_cds_seq_segments
Inputs       : Nucleotide sequence string
Outputs      : Predicted CDS sequence segments and their lengths
Description  : Predicts CDS sequence segments in both forward and reverse directions
'''
def _predict_cds_seq_segments(seq_str):
    forward_cds_seq_segments = _generate_cds_seq_segments(seq_str, False)

    forward_seq = Seq(seq_str)
    rev_comp_seq_str = str(forward_seq.reverse_complement())
    rev_comp_cds_seq_segments = _generate_cds_seq_segments(rev_comp_seq_str, True)

    candidate_cds_seq_segments = forward_cds_seq_segments + rev_comp_cds_seq_segments

    candidate_cds_lens = set()
    for cds_seq_segment in candidate_cds_seq_segments:
        candidate_cds_lens.add(cds_seq_segment.length)

    return candidate_cds_seq_segments, candidate_cds_lens

'''
Function name: _predict_cds_lens_by_protein_seq
Inputs       : Protein sequence string
Outputs      : Lengths of predicted CDS sequence segments
Description  : Predicts CDS sequence segment lengths from input protein sequence string
'''
def _predict_cds_lens_by_protein_seq(seq_str):
    orf_complete_cds_len = len(seq_str) * 3
    complete_cds_len = orf_complete_cds_len + 3

    candidate_cds_lens = set()

    '''i starts from -1 to handle the last amino acid translated from partial codon'''
    for i in range(-1, 3):
        candidate_cds_lens.add(orf_complete_cds_len + i)
        candidate_cds_lens.add(complete_cds_len + i)

    return candidate_cds_lens

'''
Function name: predict_cds_regions
Inputs       : DNA and protein sequence records with nucleotide accession numbers
Outputs      : CDS region lengths as sequence length filters, predicted CDS sequence segments (for DNA
               sequence records only)
Description  : Predicts CDS sequence segments and their lengths for DNA sequence records, and the
               lengths of the CDS sequence segments for protein sequence records
'''
def predict_cds_regions(nt_id_nt_seq_records, nt_id_protein_seq_records):
    cds_seq_len_filters = dict()
    candidate_cds_seq_segment_map = dict()

    for nt_acc_num, nt_seq_records in nt_id_nt_seq_records.items():
        nt_non_ver_acc_num = trim_version(nt_acc_num)

        for seq_record in nt_seq_records:
            candidate_cds_seq_segments, candidate_cds_lens = _predict_cds_seq_segments(str(seq_record.seq))
            add_to_group(candidate_cds_seq_segment_map, seq_record.id, candidate_cds_seq_segments)

            if nt_non_ver_acc_num in cds_seq_len_filters:
                cds_seq_len_filters[nt_non_ver_acc_num].update(candidate_cds_lens)
            else:
                cds_seq_len_filters[nt_non_ver_acc_num] = candidate_cds_lens

    for nt_acc_num, protein_seq_records in nt_id_protein_seq_records.items():
        nt_non_ver_acc_num = trim_version(nt_acc_num)

        for seq_record in protein_seq_records:
            candidate_cds_lens = _predict_cds_lens_by_protein_seq(str(seq_record.seq))
            if nt_non_ver_acc_num in cds_seq_len_filters:
                cds_seq_len_filters[nt_non_ver_acc_num].update(candidate_cds_lens)
            else:
                cds_seq_len_filters[nt_non_ver_acc_num] = candidate_cds_lens

    return cds_seq_len_filters, candidate_cds_seq_segment_map
