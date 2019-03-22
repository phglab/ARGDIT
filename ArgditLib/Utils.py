'''Module for the utility functions used by ARGDIT tools'''

from .Constants import NT_ACC_NUM_EMBED_PATTERN, PROTEIN_ACC_NUM_EMBED_PATTERN, WGS_ACC_NUM_PATTERN, \
     IUPAC_DNA_STR_PATTERN, IUPAC_PROTEIN_STR_PATTERN, IUPAC_AMBIG_DNA_BASES, IUPAC_AMBIG_PROTEIN_BASES, \
     NCBI_LIVE_SEQ_STATUS
from .ProcLog import ProcLog
from functools import partial
import os
import re

'''
Function name: check_acc_num_type
Inputs       : FASTA sequence header
Outputs      : NCBI accession number type (NT/AA/empty string)
Description  : First detects for the NCBI nucleotide accession number, and then detects for the NCBI
               protein accession number when no nucleotide accession number is present. Returns an empty
               string when no accession number is detected
'''
def check_acc_num_type(fasta_header):
    is_acc_num_found = False

    m = re.search(NT_ACC_NUM_EMBED_PATTERN, fasta_header)
    if m:
        nt_acc_num_start_pos = m.start(2)
        is_acc_num_found = True
    else:
        nt_acc_num_start_pos = len(fasta_header)

    m = re.search(PROTEIN_ACC_NUM_EMBED_PATTERN, fasta_header)
    if m:
        protein_acc_num_start_pos = m.start(2)
        is_acc_num_found = True
    else:
        protein_acc_num_start_pos = len(fasta_header)

    if is_acc_num_found:
        if nt_acc_num_start_pos < protein_acc_num_start_pos:
            return 'NT'
        elif nt_acc_num_start_pos > protein_acc_num_start_pos:
            return 'AA'

    return ''

'''
Function name: _extract_genbank_acc_num
Inputs       : FASTA sequence header and pattern for the accession number to be extracted
Outputs      : Accession number extracted or None
Description  : Extracts the accession number matching the pattern, or returns None when the pattern
               is not matched
'''
def _extract_genbank_acc_num(fasta_header, acc_num_pattern):
    m = re.search(acc_num_pattern, fasta_header)
    if m:
        return m.group(2)
    else:
        return None

'''
Function name: extract_nt_acc_num
Inputs       : FASTA sequence header
Outputs      : Nucleotide accession number or None
Description  : Invokes _extract_genbank_acc_num with nucleotide accession number pattern
'''
extract_nt_acc_num = partial(_extract_genbank_acc_num, acc_num_pattern = NT_ACC_NUM_EMBED_PATTERN)

'''
Function name: extract_protein_acc_num
Inputs       : FASTA sequence header
Outputs      : Protein accession number or None
Description  : Invokes _extract_genbank_acc_num with protein accession number pattern
'''
extract_protein_acc_num = partial(_extract_genbank_acc_num, acc_num_pattern = PROTEIN_ACC_NUM_EMBED_PATTERN)

'''
Function name: trim_version
Inputs       : Accession number
Outputs      : Accession number without version
Description  : Trims the version number (.x) from the accesion number
'''
def trim_version(acc_num):
    acc_num_fields = acc_num.split('.')

    return acc_num_fields[0]

'''
Function name: is_non_version_acc_num
Inputs       : Accession number
Outputs      : True or False
Description  : Checks if the accession number has not version
'''
def is_non_version_acc_num(acc_num):
    return (acc_num.find('.') == -1)

'''
Function name: split_wgs_acc_num
Inputs       : Nucleotide accession numbers
Outputs      : WGS accession number list and non-WGS accession number list
Description  : Extracts WGS accession numbers from the input nucleotide accession numbers
'''
def split_wgs_acc_nums(nt_acc_nums):
    wgs_acc_nums = list()
    non_wgs_acc_nums = list()

    for nt_acc_num in nt_acc_nums:
        if re.match(WGS_ACC_NUM_PATTERN, nt_acc_num):
            wgs_acc_nums.append(nt_acc_num)
        else:
            non_wgs_acc_nums.append(nt_acc_num)

    return (wgs_acc_nums, non_wgs_acc_nums)

'''
Function name: check_seq_type
Inputs       : Sequence string
Outputs      : Sequence type (NT/AA/empty string)
Description  : Checks and returns the sequence type for nucleotide (NT), protein (AA), or unknown (empty
               string). Using IUPAC codes, an input sequence is said to be a nucleotide sequence when at
               least 95% of its bases are unambigious (i.e. must be ACGT). The same applies to protein
               sequence
'''
def check_seq_type(seq_str):
    if re.match(IUPAC_DNA_STR_PATTERN, seq_str):
        unambig_dna_seq_str = re.sub(IUPAC_AMBIG_DNA_BASES, '', seq_str)
        if len(unambig_dna_seq_str) / len(seq_str) >= 0.95:
            return 'NT'
        else:
            return 'AA'
    elif re.match(IUPAC_PROTEIN_STR_PATTERN, seq_str):
        unambig_protein_seq_str = re.sub(IUPAC_AMBIG_PROTEIN_BASES, '', seq_str)
        if len(unambig_protein_seq_str) / len(seq_str) >= 0.95:
            return 'AA'

    return ''

'''
Function name: check_seq_completeness
Inputs       : Booleans indicating 5' partial or 3' partial status
Outputs      : Complete/5_partial/3_partial/5_and 3_partial
Description  : Determines the sequence completeness status according to the two input booleans
'''
def check_seq_completeness(is_5_partial, is_3_partial):
    if is_5_partial is None or is_3_partial is None:
        return ''

    if is_5_partial:
        if is_3_partial:
            return '5_and_3_partial'
        else:
            return '5_partial'
    else:
        if is_3_partial:
            return '3_partial'
        else:
            return 'complete'

'''
Function name: rectify_annotation_header_field
Inputs       : Sequence header field
Outputs      : Cleansed sequence header field
Description  : Removes some special characters from the header field, and replaces any space
               character by an underscore '_'
'''
def rectify_annotation_header_field(header_field):
    spaced_header_field = re.sub(r'[^a-zA-Z0-9_\- ]', r'', header_field)

    return re.sub(r'[ ]+', r'_', spaced_header_field)

'''
Function name: rectify_gene_short_name
Inputs       : Gene name in short form
Outputs      : Cleansed gene name
Description  : Removes some special characters and space character from the gene name
'''
def rectify_gene_short_name(gene_short_name):
    return re.sub(r'[^a-zA-Z0-9_\-]', r'', gene_short_name)

'''
Function name: create_supp_file_path
Inputs       : Original file path and the supplementary file name suffix
Outputs      : File path
Description  : Creates a new supplementary file name for the original file name by removing the file
               extension from it, and then appends with the input supplementary file name suffix. The
               file path for the created supplementary file name is returned.
'''
def create_supp_file_path(org_file_path, supp_file_suffix):
    org_file_name = os.path.basename(org_file_path)
    file_path_comps = org_file_name.rpartition('.')
    if file_path_comps[1] == '':
        supp_file_name = '{}{}'.format(org_file_name, supp_file_suffix)
    else:
        supp_file_name = '{}{}'.format(file_path_comps[0], supp_file_suffix)

    return os.path.join(os.path.dirname(org_file_path), supp_file_name)

'''
Function name: match_non_acc_fmt_genbank_protein_seqs
Inputs       : Candidate translated protein sequences (strings), protein information downloaded from
               NCBI, and protein IDs not in NCBI protein accession number format
Outputs      : Protein IDs (not in NCBI protein accession number format)
Description  : It may happen that the identifier for the downloaded protein information is not the
               protein accession number. Such identifiers are first gathered so that protein
               sequence matching can still be performed with this function
'''
def match_non_acc_fmt_genbank_protein_seqs(candidate_protein_seq_strs, genbank_protein_info_set,
                                           non_acc_fmt_protein_ids):
    if type(candidate_protein_seq_strs) is str:
        candidate_protein_seq_strs = {candidate_protein_seq_strs}

    for protein_id in non_acc_fmt_protein_ids:
        protein_info = genbank_protein_info_set[protein_id]
        if protein_info.seq_str in candidate_protein_seq_strs:
            return protein_id

    return None

'''
Function name: add_to_group
Inputs       : Existing data groups, target group ID, and input data
Outputs      : Nil
Description  : Adds input data to the target group. Create a new group when the target group
               does not exist
'''
def add_to_group(data_groups, group_id, input_data):
    input_data_type = type(input_data)

    if input_data_type is list:
        if group_id in data_groups:
            data_groups[group_id] += input_data
        else:
            data_groups[group_id] = list(input_data)
    elif input_data_type is set:
        if group_id in data_groups:
            data_groups[group_id].update(input_data)
        else:
            data_groups[group_id] = set(input_data)
    else:
        if group_id in data_groups:
            data_groups[group_id].append(input_data)
        else:
            data_groups[group_id] = [input_data]

'''
Function name: extract_seq_class_by_fields
Inputs       : Sequence header, sequence class field indices, original field separator, new
               field separator
Outputs      : Extracted sequence class or None
Description  : Extracts the sequence class fields from the sequence header, and then joins
               them to form a sequence class label, in which the fields are separated by a
               new field separator
'''
def extract_seq_class_by_fields(fasta_header, extract_index_list, org_sep, new_sep):
    header_fields = fasta_header.split(org_sep)
    header_field_range =  set(range(-1 * len(header_fields), len(header_fields)))
    if not (set(extract_index_list) <= header_field_range):
        return None

    seq_class_label = None

    for extract_index in extract_index_list:
        if seq_class_label is None:
            seq_class_label = header_fields[extract_index]
        else:
            seq_class_label = '{}{}{}'.format(seq_class_label, new_sep, header_fields[extract_index])

    return seq_class_label

'''
Function name: group_protein_by_seq_class
Inputs       : Existing protein groups, protein sequence record, sequence class field indices,
               and ARGDIT configuration
Outputs      : Nil
Description  : Categorizes protein sequence records according to their extracted sequence
               class labels
'''
def group_protein_by_seq_class(seq_class_protein_seq_record_grps, protein_seq_record, class_label_field_nums, config):
    seq_class_label = extract_seq_class_by_fields(protein_seq_record.id, class_label_field_nums,
                                                  config.fasta_header_field_sep, config.field_sep)
    if seq_class_label is None:
        ProcLog.log_exec_error('Sequence class cannot be extracted from {}'.format(protein_seq_record.description))
    else:
        add_to_group(seq_class_protein_seq_record_grps, seq_class_label, protein_seq_record)

'''
Function name: check_redundant_seq
Inputs       : Sequence, sequence ID, validated sequences (strings), and input sequence type
               (optional)
Outputs      : A 2-tuple showing an existing examined sequence ID and reverse complementarity
               flag, or a 2-tuple of None
Description  : Checks whether the input sequence is redundant or reverse complemented with
               any of the existing examined sequences
'''
def check_redundant_seq(seq, seq_id, validated_seq_strs, seq_type = None):
    forward_seq_str = str(seq)
    if seq_type is None:
        seq_type = check_seq_type(forward_seq_str)

    if seq_type == 'NT':
        rev_comp_seq_str = str(seq.reverse_complement())
    else:
        rev_comp_seq_str = ''

    if forward_seq_str in validated_seq_strs:
        return validated_seq_strs[forward_seq_str], False
    elif rev_comp_seq_str in validated_seq_strs:
        return validated_seq_strs[rev_comp_seq_str], True
    else:
        validated_seq_strs[forward_seq_str] = seq_id
        return None, None

'''
Function name: is_duplicated_item
Inputs       : Item to check and existing examined items
Outputs      : True or False
Description  : Checks whether the target item already exists in the existing examined items,
               and adds it to the examined items when absent
'''
def is_duplicated_item(item_to_check, checked_items):
    if item_to_check in checked_items:
        return True
    else:
        checked_items.add(item_to_check)
        return False

'''
Function name: is_obsolete_ncbi_seq
Input        : sequence status
Output       : True or False
Description  : Checks whether the input sequence status refers to a live sequence in the
               NCBI repository
'''
def is_obsolete_ncbi_seq(seq_status):
    return seq_status != NCBI_LIVE_SEQ_STATUS
