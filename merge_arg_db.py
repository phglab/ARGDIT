#!/usr/bin/env python3

from ArgditLib import CDSPredict
from ArgditLib import EntrezDBAccess
from ArgditLib import OptionParser
from ArgditLib import Translate
from ArgditLib import Utils
from ArgditLib.Config import Config
from ArgditLib.Constants import SEQ_VERSION_MARKUP
from ArgditLib.OntologyAnnotator import OntologyAnnotator
from ArgditLib.ProcLog import ProcLog
from ArgditLib.RepositoryBuffer import RepositoryBuffer
from ArgditLib.SequenceFileParser import SequenceFileParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
import sys

'''ARG database integration main program'''

'''
Function name: match_protein_cds_for_nt_id_nt_seq_rec
Inputs       : Nucleotide sequences with nucleotide accession numbers, ARG database file path,
               downloaded NCBI data, nucleotide accession number information, boolean indicating
               whether the ARG database is a schema database
Outputs      : Translated protein sequences for the nucleotide sequences in the ARG database, CDS
               regions facilitating correct translation for the nucleotide sequences
Description  : Match the ARG nucleotide sequences by translating them with the associated CDS regions
               (provided by the downloaded NCBI data), and then matching the translation products with
               the protein sequences referred by the CDS regions
'''
def match_protein_cds_for_nt_id_nt_seq_rec(nt_id_nt_seq_records, seq_db_path, data_params, nt_acc_num_params,
                                           is_schema_db):
    (genbank_protein_info_set, candidate_cds_seq_segment_map, target_cds_region_grps,
     non_acc_fmt_protein_ids) = data_params
    latest_ver_nt_acc_num_map, matched_nt_non_ver_acc_nums = nt_acc_num_params
    
    matched_protein_info_set = dict()
    matched_cds_region_set = dict()

    for nt_acc_num, nt_seq_records in nt_id_nt_seq_records.items():
        nt_non_ver_acc_num = Utils.trim_version(nt_acc_num)

        if nt_non_ver_acc_num not in matched_nt_non_ver_acc_nums:
            for nt_seq_record in nt_seq_records:
                ProcLog.log_acc_num_not_found(msg = nt_seq_record.id, db_name = seq_db_path,
                                              is_for_schema_db = is_schema_db)

            continue

        latest_ver_nt_acc_num = latest_ver_nt_acc_num_map[nt_non_ver_acc_num]

        for nt_seq_record in nt_seq_records:
            '''
            Try to translate the predicted CDS sequences using the target CDS regions, and then compare
            with the protein products specified in the target CDS regions. Matched protein is stored in
            translated_protein_info and matched_cds_seq_region contains its associated CDS region
            '''
            translated_protein_info, _, matched_target_cds_region, _ = \
                Translate.search_correct_cds_translation(candidate_cds_seq_segment_map[nt_seq_record.id],
                                                         target_cds_region_grps[latest_ver_nt_acc_num],
                                                         genbank_protein_info_set, non_acc_fmt_protein_ids)
            if translated_protein_info is None:
                ProcLog.log_seq_mismatch(seq_id = nt_seq_record.id, is_obsolete_ver = False, db_name = seq_db_path,
                                         is_for_schema_db = is_schema_db)
            else:
                matched_protein_info_set[nt_seq_record.id] = translated_protein_info
                matched_cds_region_set[nt_seq_record.id] = matched_target_cds_region

    return matched_protein_info_set, matched_cds_region_set

'''
Function name: match_protein_cds_for_nt_id_protein_seq_rec
Inputs       : Protein sequences with nucleotide accession numbers, ARG database file path, downloaded
               NCBI data, nucleotide accession number information, boolean indicating whether the ARG
               database is a schema database
Outputs      : Matched ARG protein sequences, CDS regions associated with the matched ARG protein
               sequences
Description  : Match the ARG protein sequences by matching them with the protein products referred by
               the associated CDS regions (provided by the downloaded NCBI data)
'''
def match_protein_cds_for_nt_id_protein_seq_rec(nt_id_protein_seq_records, seq_db_path, data_params,
                                                nt_acc_num_params, is_schema_db):
    genbank_protein_info_set, _, target_cds_region_grps, non_acc_fmt_protein_ids = data_params
    latest_ver_nt_acc_num_map, matched_nt_non_ver_acc_nums = nt_acc_num_params
    
    matched_protein_info_set = dict()
    matched_cds_region_set = dict()

    for nt_acc_num, protein_seq_records in nt_id_protein_seq_records.items():
        nt_non_ver_acc_num = Utils.trim_version(nt_acc_num)

        if nt_non_ver_acc_num not in matched_nt_non_ver_acc_nums:
            for protein_seq_record in protein_seq_records:
                ProcLog.log_acc_num_not_found(msg = protein_seq_record.id, db_name = seq_db_path,
                                              is_for_schema_db = is_schema_db)

            continue

        latest_ver_nt_acc_num = latest_ver_nt_acc_num_map[nt_non_ver_acc_num]

        for protein_seq_record in protein_seq_records:
            protein_seq_str = str(protein_seq_record.seq)
            matched_target_cds_region = None

            '''
            Compare the ARG protein sequence with the protein products specified in the target CDS
            regions
            '''
            for target_cds_region in target_cds_region_grps[latest_ver_nt_acc_num]:
                if target_cds_region.protein_id in genbank_protein_info_set:
                    genbank_protein_info = genbank_protein_info_set[target_cds_region.protein_id]
                    if protein_seq_str == genbank_protein_info.seq_str:
                        matched_target_cds_region = target_cds_region
                        break
                else:
                    '''
                    Search the ARG protein sequence in protein products identified by non-accession
                    format identifiers
                    '''
                    matched_protein_id = Utils.match_non_acc_fmt_genbank_protein_seqs(protein_seq_str,
                                                                                      genbank_protein_info_set,
                                                                                      non_acc_fmt_protein_ids)
                    if matched_protein_id is not None:
                        genbank_protein_info = genbank_protein_info_set[matched_protein_id]
                        matched_target_cds_region = target_cds_region
                        break

            if matched_target_cds_region is not None:
                matched_protein_info_set[protein_seq_record.id] = genbank_protein_info
                matched_cds_region_set[protein_seq_record.id] = matched_target_cds_region
            else:
                ProcLog.log_seq_mismatch(seq_id = protein_seq_record.id, is_obsolete_ver = False, db_name = seq_db_path,
                                         is_for_schema_db = is_schema_db)

    return matched_protein_info_set, matched_cds_region_set

'''
Function name: match_protein_info_for_protein_id_nt_seq_rec
Inputs       : Nucleotide sequences with protein accession numbers, ARG database file path, downloaded
               NCBI data, protein accession number information, boolean indicating whether the ARG
               database is a schema database
Outputs      : Translated protein sequences for the nucleotide sequences in the ARG database
Description  : Match the ARG nucleotide sequences by performing six-frame translation and then matching
               the translation products for each frame with the protein sequences referred by the
               protein accession numbers
'''
def match_protein_info_for_protein_id_nt_seq_rec(protein_id_nt_seq_records, seq_db_path, data_params,
                                                 protein_acc_num_params, is_schema_db):
    genbank_protein_info_set, _, _, non_acc_fmt_protein_ids = data_params
    latest_ver_protein_acc_num_map, matched_protein_non_ver_acc_nums = protein_acc_num_params

    matched_protein_info_set = dict()

    for protein_acc_num, nt_seq_records in protein_id_nt_seq_records.items():
        protein_non_ver_acc_num = Utils.trim_version(protein_acc_num)

        if protein_non_ver_acc_num in matched_protein_non_ver_acc_nums:
            latest_ver_protein_acc_num = latest_ver_protein_acc_num_map[protein_non_ver_acc_num]
        else:
            latest_ver_protein_acc_num = None
       
        for nt_seq_record in nt_seq_records:
            '''Translate the ARG nucleotide sequence and obtain products from all six frames'''
            translated_protein_seq_strs = Translate.translate(str(nt_seq_record.seq))

            if latest_ver_protein_acc_num is not None:
                genbank_protein_info = genbank_protein_info_set[latest_ver_protein_acc_num]
                is_protein_seq_matched = genbank_protein_info.seq_str in translated_protein_seq_strs
            else:
                '''
                Search the translation candidates in protein products identified by non-accession
                format identifiers
                '''
                matched_protein_id = Utils.match_non_acc_fmt_genbank_protein_seqs(translated_protein_seq_strs,
                                                                                  genbank_protein_info_set,
                                                                                  non_acc_fmt_protein_ids)
                if matched_protein_id is None:
                    is_protein_seq_matched = False
                else:
                    genbank_protein_info = genbank_protein_info_set[matched_protein_id]
                    is_protein_seq_matched = True

            if is_protein_seq_matched:
                matched_protein_info_set[nt_seq_record.id] = genbank_protein_info
            else:
                ProcLog.log_seq_mismatch(seq_id = nt_seq_record.id, is_obsolete_ver = False, db_name = seq_db_path,
                                         is_for_schema_db = is_schema_db)

    return matched_protein_info_set

'''
Function name: match_protein_info_for_protein_id_protein_seq_rec
Inputs       : Protein sequences with protein accession numbers, ARG database file path, downloaded
               NCBI data, protein accession number information, boolean indicating whether the ARG
               database is a schema database
Outputs      : Matched ARG protein sequences
Description  : Match the ARG protein sequences by matching them with the protein sequences referred by
               the protein accession numbers
'''
def match_protein_info_for_protein_id_protein_seq_rec(protein_id_protein_seq_records, seq_db_path, data_params,
                                                      protein_acc_num_params, is_schema_db):
    genbank_protein_info_set, _, _, non_acc_fmt_protein_ids = data_params
    latest_ver_protein_acc_num_map, matched_protein_non_ver_acc_nums = protein_acc_num_params

    matched_protein_info_set = dict()

    for protein_acc_num, protein_seq_records in protein_id_protein_seq_records.items():
        protein_non_ver_acc_num = Utils.trim_version(protein_acc_num)

        if protein_non_ver_acc_num in matched_protein_non_ver_acc_nums:
            latest_ver_protein_acc_num = latest_ver_protein_acc_num_map[protein_non_ver_acc_num]
        else:
            latest_ver_protein_acc_num = None

        for protein_seq_record in protein_seq_records:
            protein_seq_str = str(protein_seq_record.seq)
        
            if latest_ver_protein_acc_num is not None:
                genbank_protein_info = genbank_protein_info_set[latest_ver_protein_acc_num]
                is_protein_seq_matched = (protein_seq_str == genbank_protein_info.seq_str)
            else:
                '''
                Search the ARG protein sequence in protein products identified by non-accession format
                identifiers
                '''
                matched_protein_id = Utils.match_non_acc_fmt_genbank_protein_seqs(protein_seq_str,
                                                                                  genbank_protein_info_set,
                                                                                  non_acc_fmt_protein_ids)
                if matched_protein_id is None:
                    is_protein_seq_matched = False
                else:
                    genbank_protein_info = genbank_protein_info_set[matched_protein_id]
                    is_protein_seq_matched = True

            if is_protein_seq_matched:
                matched_protein_info_set[protein_seq_record.id] = genbank_protein_info
            else:
                ProcLog.log_seq_mismatch(seq_id = protein_seq_record.id, is_obsolete_ver = False, db_name = seq_db_path,
                                         is_for_schema_db = is_schema_db)

    return matched_protein_info_set

'''
Function name: match_protein_cds
Inputs       : Sequence file parser, ARG database file path, downloaded NCBI data, nucleotide accession
               number information, protein accession number information, boolean indicating whether the
               ARG database is a schema database, boolean controlling whether to process nucleotide
               sequences only
Outputs      : Translated protein sequences for the nucleotide sequences in the ARG database, CDS
               regions facilitating correct translation for the nucleotide sequences
Description  : Match the ARG nucleotides/protein sequences identified by either nucleotide or protein
               accession numbers, call different validation functions according to each of the four
               sequence-identifier categories (i.e. nucleotide sequences with nucleotide accession
               numbers, nucleotide sequences with protein accession numbers, etc.)
'''
def match_protein_cds(seq_file_parser, seq_db_path, data_params, nt_acc_num_params, protein_acc_num_params,
                      is_schema_db = False, is_nt_seq_only = False):
    all_matched_protein_info_set = dict()
    all_matched_cds_region_set = dict()

    '''Nucleotide sequences with nucleotide accession numbers'''
    nt_id_nt_seq_records = seq_file_parser.get_nt_id_nt_seq_records()
    if len(nt_id_nt_seq_records) > 0:
        matched_protein_info_set, matched_cds_region_set = \
            match_protein_cds_for_nt_id_nt_seq_rec(nt_id_nt_seq_records, seq_db_path, data_params, nt_acc_num_params,
                                                   is_schema_db)
        all_matched_protein_info_set.update(matched_protein_info_set)
        all_matched_cds_region_set.update(matched_cds_region_set)

    '''Protein sequences with nucleotide accession numbers'''
    nt_id_protein_seq_records = seq_file_parser.get_nt_id_protein_seq_records()
    if not is_nt_seq_only and len(nt_id_protein_seq_records) > 0:
        matched_protein_info_set, matched_cds_region_set = \
            match_protein_cds_for_nt_id_protein_seq_rec(nt_id_protein_seq_records, seq_db_path, data_params,
                                                        nt_acc_num_params, is_schema_db)
        all_matched_protein_info_set.update(matched_protein_info_set)
        all_matched_cds_region_set.update(matched_cds_region_set)

    '''Nucleotide sequences with protein accession numbers'''
    protein_id_nt_seq_records = seq_file_parser.get_protein_id_nt_seq_records()
    if len(protein_id_nt_seq_records) > 0:
        matched_protein_info_set = \
            match_protein_info_for_protein_id_nt_seq_rec(protein_id_nt_seq_records, seq_db_path, data_params,
                                                         protein_acc_num_params, is_schema_db)
        all_matched_protein_info_set.update(matched_protein_info_set)

    '''Protein sequences with protein accession numbers'''
    protein_id_protein_seq_records = seq_file_parser.get_protein_id_protein_seq_records()
    if not is_nt_seq_only and len(protein_id_protein_seq_records) > 0:
        matched_protein_info_set = \
            match_protein_info_for_protein_id_protein_seq_rec(protein_id_protein_seq_records, seq_db_path,
                                                              data_params, protein_acc_num_params, is_schema_db)
        all_matched_protein_info_set.update(matched_protein_info_set)

    return all_matched_protein_info_set, all_matched_cds_region_set

'''
Function name: build_ontology_annotator
Inputs       : Sequence file parser, translated protein sequences for the nucleotide sequences in the
               ARG schema database, ARG ontology label field numbers, ARGDIT configuration
Outputs      : Ontology annotator
Description  : Generate the ontology annotator using the ARG schema database which is parsed by the
               corresponding sequence file parser. The ARG ontology class information in the sequence
               header is indicated by the ontology label field numbers
'''
def build_ontology_annotator(seq_file_parser, matched_protein_info_set, otl_label_field_nums, config):
    otl_protein_seq_record_grps = dict()

    for nt_seq_record in seq_file_parser.get_nt_seq_record_list():
        if nt_seq_record.id not in matched_protein_info_set:
            continue

        matched_protein_info = matched_protein_info_set[nt_seq_record.id]
        protein_seq_record = SeqRecord(Seq(matched_protein_info.seq_str), id = nt_seq_record.id, name = '',
                                       description = '')
        Utils.group_protein_by_otl_class(otl_protein_seq_record_grps, protein_seq_record, otl_label_field_nums,
                                         config)

    for protein_seq_record in seq_file_parser.get_protein_seq_record_list():
        Utils.group_protein_by_otl_class(otl_protein_seq_record_grps, protein_seq_record, otl_label_field_nums,
                                         config)

    '''Ontology annotator only includes ARG ontology class with more than config.min_seq_count sequences'''
    qualified_otl_seq_record_grps = dict()

    for otl_class, seq_records in otl_protein_seq_record_grps.items():
        if len(seq_records) >= config.min_seq_count:
            qualified_otl_seq_record_grps[otl_class] = seq_records

    return OntologyAnnotator(qualified_otl_seq_record_grps), set(otl_protein_seq_record_grps.keys())

'''
Function name: batch_otl_annotate
Inputs       : Ontology annotator, nucleotide sequences to be annotated, protein sequences to be
               annotated, translated protein sequences for the nucleotide sequences in the ARG
               databases
Outputs      : ARG ontology class labels
Description  : Predict ARG ontology class labels for the input nucleotide and/or protein sequences
               using the ontology annotator
'''
def batch_otl_annotate(otl_annotator, nt_seq_record_list, protein_seq_record_list, matched_protein_info_set):
    input_protein_seq_record_list = list(protein_seq_record_list)
    for nt_seq_record in nt_seq_record_list:
        if nt_seq_record.id in matched_protein_info_set:
            matched_protein_info = matched_protein_info_set[nt_seq_record.id]
            input_protein_seq_record_list.append(SeqRecord(Seq(matched_protein_info.seq_str),
                                                           id = nt_seq_record.id, name = nt_seq_record.name,
                                                           description = nt_seq_record.description))

    return otl_annotator.annotate_ontology_label(input_protein_seq_record_list)

'''
Function name: generate_genbank_seq_header
Inputs       : ARG configuration, accession number, information for the protein sequences downloaded
               from the NCBI protein database, (if any) associated CDS region for the ARG nucleotide
               sequence
Outputs      : Sequence header
Description  : Generate sequence header using the download protein information and (if any) associated
               CDS region information
'''
def generate_genbank_seq_header(config, acc_num, matched_protein_info, matched_target_cds_region = None):
    protein_name = Utils.rectify_annotation_header_field(matched_protein_info.name)
    organism_name = Utils.rectify_annotation_header_field(matched_protein_info.organism)
    field_sep = config.fasta_header_field_sep

    if matched_target_cds_region is None:
        if matched_protein_info.coding_gene_short_name is None:
            coding_gene_short_name = ''
        else:
            coding_gene_short_name = Utils.rectify_gene_short_name(matched_protein_info.coding_gene_short_name)

        seq_header_template = '{}' * 11
        
        return seq_header_template.format(acc_num, field_sep, coding_gene_short_name, field_sep, protein_name,
                                          field_sep, SEQ_VERSION_MARKUP, field_sep, organism_name,
                                          field_sep * 3, matched_protein_info.seq_completeness())
    else:
        if matched_target_cds_region.gene_short_name is None:
            gene_short_name = ''
        else:
            gene_short_name = Utils.rectify_gene_short_name(matched_target_cds_region.gene_short_name)

        seq_header_template = '{}' * 15

        return seq_header_template.format(acc_num, field_sep, gene_short_name, field_sep, protein_name,
                                          field_sep, SEQ_VERSION_MARKUP, field_sep, organism_name,
                                          field_sep, matched_target_cds_region.get_region_range_annotation(),
                                          field_sep, matched_target_cds_region.length, field_sep,
                                          matched_target_cds_region.seq_completeness())

'''Entry point of the main program'''
parser = argparse.ArgumentParser()
parser.add_argument('seq_db_paths', nargs = '+', help = 'nucleotide/protein database FASTA file paths')
parser.add_argument('-o', '--output', action = 'store', dest = 'output_seq_db_path', required = True,
                    help = 'output database file path')
parser.add_argument('-s', '--schema', action = 'store', dest = 'otl_schema', nargs = 2,
                    help = 'schema database and ontology field numbers for sequence classification')
parser.add_argument('-a', '--annotate', action = 'store_true',
                    help = 'automatic re-annotation using NCBI database information')
parser.add_argument('-p', '--protein', action = 'store_true', help = 'export protein sequences')
parser.add_argument('-r', '--redundant', action = 'store_true', help = 'allow redundant sequences')
parser.add_argument('-e', '--exportlog', action = 'store_true', help = 'export integration results and process log')

args = parser.parse_args()

ProcLog.init_logs()

config = Config('config.ini')

for seq_db_path in args.seq_db_paths:
    if not os.path.exists(seq_db_path):
        ProcLog.log_exec_error('Database file \'{}\' does not exist'.format(seq_db_path))

otl_schema_file_path = None
if args.otl_schema is not None:
    otl_label_field_nums = OptionParser.parse_ontology_label_field_nums(args.otl_schema[1])

    if os.path.exists(args.otl_schema[0]):
        otl_schema_file_path = args.otl_schema[0]
    else:
        ProcLog.log_exec_error('Schema database file \'{}\' does not exist'.format(args.otl_schema[0]))

if ProcLog.has_exec_error():
    ProcLog.export_exec_error(sys.stdout)
    sys.exit()

db_seq_file_parsers = dict()
source_db_seq_rec_count = 0

for seq_db_path in args.seq_db_paths:
    if seq_db_path not in db_seq_file_parsers:
        seq_file_parser = SequenceFileParser()
        seq_file_parser.parse(seq_db_path, is_ignore_acc_num_type = not args.annotate)
        db_seq_file_parsers[seq_db_path] = seq_file_parser
        source_db_seq_rec_count += seq_file_parser.get_seq_record_count()

schema_db_seq_rec_count = 0

if args.otl_schema is not None:
    if otl_schema_file_path in db_seq_file_parsers:
        otl_schema_seq_file_parser = db_seq_file_parsers[otl_schema_file_path]
    else:
        '''
        If the ARG schema database is not one of the source ARG databases to be merged, then only its
        nucleotide sequences require NCBI database query to perform translation, and so the sequence
        file parser only parses these sequences
        '''
        otl_schema_seq_file_parser = SequenceFileParser()
        otl_schema_seq_file_parser.parse_ontology_schema_db(otl_schema_file_path, is_nt_seq_only = True)
        db_seq_file_parsers[otl_schema_file_path] = otl_schema_seq_file_parser

    schema_db_seq_rec_count += otl_schema_seq_file_parser.get_seq_record_count()

global_cds_seq_len_filters = dict()
db_class_candidate_cds_seq_segments = dict()
query_protein_acc_nums = set()

'''
NCBI database query required only when auto-generation of ARG sequence header, conversion to protein
sequences, or ARG ontology annotation is required
'''
if args.annotate or args.protein or args.otl_schema is not None:
    '''
    For ARG nucleotide and protein sequences with NCBI nucleotide accession numbers, predict the
    lengths of their potential CDS sequences, and use these lengths as sequence length filters to
    select target CDS regions; the predicted CDS sequences are also stored for ARG nucleotide
    sequences
    '''
    for seq_db_path, seq_file_parser in db_seq_file_parsers.items():
        nt_id_nt_seq_records = seq_file_parser.get_nt_id_nt_seq_records()

        if args.annotate or args.otl_schema is not None:
            nt_id_protein_seq_records = seq_file_parser.get_nt_id_protein_seq_records()
        else:
            nt_id_protein_seq_records = dict()

        if len(nt_id_nt_seq_records) > 0 or len(nt_id_protein_seq_records) > 0:
            cds_seq_len_filters, candidate_cds_seq_segment_map = \
                CDSPredict.predict_cds_regions(nt_id_nt_seq_records, nt_id_protein_seq_records)
            for nt_non_ver_acc_num, candidate_cds_lens in cds_seq_len_filters.items():
                Utils.add_to_group(global_cds_seq_len_filters, nt_non_ver_acc_num, candidate_cds_lens)

            db_class_candidate_cds_seq_segments[seq_db_path] = candidate_cds_seq_segment_map

        '''Add the protein accession numbers of the ARG sequences to the NCBI protein query'''
        protein_id_nt_seq_records = seq_file_parser.get_protein_id_nt_seq_records()
        if len(protein_id_nt_seq_records) > 0:
            query_protein_acc_nums.update(map(Utils.trim_version, protein_id_nt_seq_records.keys()))

        if args.annotate or args.otl_schema is not None:
            protein_id_protein_seq_records = seq_file_parser.get_protein_id_protein_seq_records()
            if len(protein_id_protein_seq_records) > 0:
                query_protein_acc_nums.update(map(Utils.trim_version, protein_id_protein_seq_records.keys()))

EntrezDBAccess.set_entrez_email(config.entrez_email)

if len(global_cds_seq_len_filters) > 0:
    print('Retrieving information from NCBI nucleotide database...')
    '''
    Target CDS regions contain potential CDS sequence matches, i.e. an ARG nucleotide/protein sequence
    should come from one of the relevant target CDS regions
    '''
    target_cds_region_grps, target_cds_protein_acc_nums, is_parse_complete = \
        EntrezDBAccess.search_target_cds_by_nt_acc_num(global_cds_seq_len_filters.keys(),
                                                       global_cds_seq_len_filters)

    if not is_parse_complete:
        ProcLog.log_data_retrieval_error()

    if ProcLog.has_exec_error():
        ProcLog.export_exec_error(sys.stdout)
        sys.exit()

    '''Add the protein accession numbers of the target CDS regions to the NCBI protein query'''
    query_protein_acc_nums.update(target_cds_protein_acc_nums)

    '''Keep the latest nucleotide accession version to determine potential obsolete sequences'''
    latest_ver_nt_acc_num_map = dict()
    for nt_acc_num in target_cds_region_grps.keys():
        latest_ver_nt_acc_num_map[Utils.trim_version(nt_acc_num)] = nt_acc_num

    '''Store (non-versioned) nucleotide accession numbers of the target CDS regions'''
    matched_nt_non_ver_acc_nums = latest_ver_nt_acc_num_map.keys()
else:
    target_cds_region_grps = dict()
    latest_ver_nt_acc_num_map = dict()
    matched_nt_non_ver_acc_nums = set()

'''
Nucleotide accession number information includes latest nucleotide accession version and non-versioned
accession numbers of the target CDS regions for the nucleotide sequences
'''    
nt_acc_num_params = latest_ver_nt_acc_num_map, matched_nt_non_ver_acc_nums

if len(query_protein_acc_nums) > 0:
    print('Retrieving information from NCBI protein database...')
    genbank_protein_info_set = EntrezDBAccess.search_protein_info(query_protein_acc_nums)
    
    if ProcLog.has_exec_error():
        ProcLog.export_exec_error(sys.stdout)
        sys.exit()

    '''
    It is possible that some retrieved records are identifed by non-accession format identifiers, hence
    they need to be extracted for special matching
    '''
    non_acc_fmt_protein_ids = set(genbank_protein_info_set.keys()) - query_protein_acc_nums

    '''Keep the latest protein accession version to determine potential obsolete sequences'''
    latest_ver_protein_acc_num_map = dict()
    for protein_acc_num in genbank_protein_info_set.keys():
        latest_ver_protein_acc_num_map[Utils.trim_version(protein_acc_num)] = protein_acc_num

    '''Store (non-versioned) protein accession numbers of the target CDS regions'''
    matched_protein_non_ver_acc_nums = latest_ver_protein_acc_num_map.keys()
else:
    genbank_protein_info_set = dict()
    non_acc_fmt_protein_ids = set()
    latest_ver_protein_acc_num_map = dict()
    matched_protein_non_ver_acc_nums = set()

'''
Protein accession number information includes latest protein accession version and non-versioned
accession numbers of the target protein sequences
'''
protein_acc_num_params = latest_ver_protein_acc_num_map, matched_protein_non_ver_acc_nums

'''Build the ontology annotator for the ARG schema database'''
if args.otl_schema is not None:
    if otl_schema_file_path in db_class_candidate_cds_seq_segments:
        candidate_cds_seq_segment_map = db_class_candidate_cds_seq_segments[otl_schema_file_path]
    else:
        candidate_cds_seq_segment_map = dict()

    '''
    Downloaded NCBI data includes protein sequences and their basic information (protein name, etc.),
    candidate CDS sequences, target CDS regions, and non-accession format identifiers for protein
    sequences
    '''
    data_params = (genbank_protein_info_set, candidate_cds_seq_segment_map, target_cds_region_grps,
                   non_acc_fmt_protein_ids)

    '''
    Match the protein sequences (previously downloaded from the NCBI database) for the nucleotide
    sequences in the ARG schema database to build the ontology annotator
    '''
    matched_protein_info_set, _ = match_protein_cds(otl_schema_seq_file_parser, otl_schema_file_path, data_params,
                                                    nt_acc_num_params, protein_acc_num_params, is_schema_db = True,
                                                    is_nt_seq_only = True)

    '''
    If the ontology schema file path is not in the sequence DB file paths, parse the schema file again
    to recover the protein sequence records, which are omitted in the first parse to avoid unnecessary
    GenBank query
    '''
    if otl_schema_file_path not in args.seq_db_paths:
        otl_schema_seq_file_parser.parse_ontology_schema_db(otl_schema_file_path)
        for seq_header in otl_schema_seq_file_parser.get_invalid_acc_num_fmt_seq_rec_ids():
            ProcLog.log_invalid_acc_num_fmt(msg = seq_header, db_name = otl_schema_file_path, is_for_schema_db = True)

        for seq_header in otl_schema_seq_file_parser.get_unknown_seq_type_seq_rec_ids():
            ProcLog.log_unknown_seq_type(msg = seq_header, db_name = otl_schema_file_path, is_for_schema_db = True)

    otl_annotator, schema_otl_classes = build_ontology_annotator(otl_schema_seq_file_parser, matched_protein_info_set,
                                                                 otl_label_field_nums, config)

repos_buffer = RepositoryBuffer()
merged_db_seq_headers = set()

'''Process the ARG sequences in each source ARG database for export'''
for seq_db_path in args.seq_db_paths:
    seq_file_parser = db_seq_file_parsers[seq_db_path]

    for seq_header in seq_file_parser.get_unknown_seq_type_seq_rec_ids():
        ProcLog.log_unknown_seq_type(msg = seq_header, db_name = seq_db_path)

    '''
    Match the protein sequences (previously downloaded from the NCBI database) when auto-generation of
    ARG sequence header, conversion to protein sequences, or ARG ontology annotation is required
    '''
    if args.annotate or args.protein or args.otl_schema is not None:
        if seq_db_path in db_class_candidate_cds_seq_segments:
            candidate_cds_seq_segment_map = db_class_candidate_cds_seq_segments[seq_db_path]
        else:
            candidate_cds_seq_segment_map = dict()

        '''
        Downloaded NCBI data includes protein sequences and their basic information (protein name,
        etc.), candidate CDS sequences, target CDS regions, and non-accession format identifiers for
        protein sequences
        '''
        data_params = (genbank_protein_info_set, candidate_cds_seq_segment_map, target_cds_region_grps,
                       non_acc_fmt_protein_ids)
        '''
        Match the protein sequences (previously downloaded from the NCBI database) for the nucleotide
        sequences in the source ARG database for auto-generation of sequence headers, protein
        sequence conversion, or ARG ontology class annotation (except the schema ARG database itself);
        associated CDS regions also matched for sequence header auto-generation
        '''
        if args.annotate:
            matched_protein_info_set, matched_cds_region_set = \
                match_protein_cds(seq_file_parser, seq_db_path, data_params, nt_acc_num_params,
                                  protein_acc_num_params)
        elif args.protein or (args.otl_schema is not None and seq_db_path != otl_schema_file_path):
            matched_protein_info_set, _ = match_protein_cds(seq_file_parser, seq_db_path, data_params,
                                                            nt_acc_num_params, protein_acc_num_params,
                                                            is_nt_seq_only = True)

        if args.otl_schema is not None and seq_db_path != otl_schema_file_path:
            otl_class_label_map = batch_otl_annotate(otl_annotator, seq_file_parser.get_nt_seq_record_list(),
                                                     seq_file_parser.get_protein_seq_record_list(),
                                                     matched_protein_info_set)
        else:
            otl_class_label_map = None

        for seq_header in seq_file_parser.get_invalid_acc_num_fmt_seq_rec_ids():
            ProcLog.log_invalid_acc_num(msg = seq_header, db_name = seq_db_path)

    seq_record_groups = seq_file_parser.get_hybrid_seq_records()

    for seq_records in seq_record_groups.values():
        for seq_record in seq_records:
            '''Skip sequences without matched protein sequences and information'''
            if (args.annotate or args.otl_schema is not None) and seq_record.id not in matched_protein_info_set:
                continue

            '''Perform sequence header auto-generation'''
            if args.annotate:
                matched_protein_info = matched_protein_info_set[seq_record.id]

                if Utils.check_acc_num_type(seq_record.id) == 'NT':
                    nt_acc_num = Utils.extract_nt_acc_num(seq_record.id)
                    matched_target_cds_region = matched_cds_region_set[seq_record.id]
                    output_seq_header = generate_genbank_seq_header(config, nt_acc_num, matched_protein_info,
                                                                    matched_target_cds_region)
                else:
                    protein_acc_num = Utils.extract_protein_acc_num(seq_record.id)
                    output_seq_header = generate_genbank_seq_header(config, protein_acc_num, matched_protein_info)
            else:
                output_seq_header = seq_record.description

            '''Convert nucleotide sequence to protein sequence for export'''
            if args.protein and Utils.check_seq_type(str(seq_record.seq)) == 'NT':
                if seq_record.id not in matched_protein_info_set:
                    continue

                matched_protein_info = matched_protein_info_set[seq_record.id]
                export_seq = Seq(matched_protein_info.seq_str)
            else:
                export_seq = seq_record.seq

            '''Perform ARG ontology classification'''
            if args.otl_schema is not None:
                org_otl_class_label = Utils.extract_ontology_class_by_fields(seq_record.id, otl_label_field_nums,
                                                                             config.fasta_header_field_sep,
                                                                             config.field_sep)
                otl_class_label = None
                if org_otl_class_label is not None and org_otl_class_label in schema_otl_classes:
                    if args.annotate:
                        otl_class_label = org_otl_class_label
                else:
                    if seq_record.id in otl_class_label_map:
                        otl_class_label = otl_class_label_map[seq_record.id].replace(config.field_sep,
                                                                                     config.fasta_header_field_sep)
                    else:
                        otl_class_label = 'UNKNOWN_ONTOLOGY_CLASS'

                if otl_class_label is not None:
                    output_seq_header = '{}{}{}'.format(output_seq_header, config.fasta_header_field_sep,
                                                        otl_class_label)

            repos_buffer.add(SeqRecord(export_seq, id = output_seq_header, name = '', description = ''),
                             seq_db_path, seq_record.id)

export_seq_rec_count = repos_buffer.export(args.output_seq_db_path, not args.redundant)

if args.exportlog:
    log_file_path = Utils.create_supp_file_path(args.output_seq_db_path, '.log')
    output_stream = open(log_file_path, 'w')
else:
    output_stream = sys.stdout

ProcLog.export_merge_db_check_logs(output_stream)

if export_seq_rec_count > 0:
    output_db_stmt = ProcLog.create_summary_stmt(export_seq_rec_count, 'exported to {}'.format(args.output_seq_db_path))
else:
    output_db_stmt = ProcLog.create_summary_stmt(export_seq_rec_count, 'exported')

export_seq_summary = [output_db_stmt]
ProcLog.export_merge_db_summary(output_stream, source_db_seq_rec_count, schema_db_seq_rec_count, export_seq_summary)

if args.exportlog:
    output_stream.close()
