from .Constants import SEQ_VERSION_MARKUP
from functools import partialmethod
import os

'''
Process log class to log and display/export process messages or errors
'''

'''Log IDs for various logs'''
INVALID_ACC_NUM_FMT = 1
ACC_NUM_NOT_FOUND = 2
'''CDS_NOT_FOUND = 3'''
FALSE_SEQ_CLASS = 4
UNKNOWN_SEQ_TYPE = 5
'''MIXED_SEQ_TYPE = 6'''
SEQ_MISMATCH = 7
DUPLICATED_HEADER = 8
REDUNDANT_SEQ = 9
OBSOLETE_VER = 10
EXEC_MSG = 11
EXEC_ERROR = 12
DB_MERGE_RESULT = 13

DEFAULT_DB_NAME = 'default_db'
SCHEMA_DB_TAG = '<schema>'

class ProcLog:
    PROTEIN_ID_MAPPING_MSG_TEMPLATE = '{}: Protein accession no. \'{}\' is possibly mapped to \'{}\''
    MERGE_OUTPUT_SEQ_SRC_MSG_TEMPLATE = 'Sequence \'{}\' exists in the following source database(s):{}'
    SEQ_SRC_DB_MSG_TEMPLATE = '\t - {} (original sequence header: {}){}'
    CROSS_DB_DUPLICATED_SEQ_HDR_MSG_TEMPLATE = 'Skipped duplicated sequence header {}'
    CROSS_DB_DUPLICATED_AUTOGEN_HDR_MSG_TEMPLATE = \
        'Skipped duplicated auto-generated sequence header {} (original header: {})'
    DATA_RETRIEVAL_ERROR = 'Data retrieval from NCBI cannot be completed. Please try again later'

    '''Initializes all logs'''
    @classmethod
    def init_logs(cls):
        cls._logs = dict()
        cls._logs[INVALID_ACC_NUM_FMT] = dict()
        cls._logs[ACC_NUM_NOT_FOUND] = dict()
        '''cls._logs[CDS_NOT_FOUND] = list()'''
        cls._logs[FALSE_SEQ_CLASS] = dict()
        cls._logs[UNKNOWN_SEQ_TYPE] = dict()
        '''cls._logs[MIXED_SEQ_TYPE] = list()'''
        cls._logs[SEQ_MISMATCH] = dict()
        cls._logs[DUPLICATED_HEADER] = dict()
        cls._logs[REDUNDANT_SEQ] = dict()
        cls._logs[OBSOLETE_VER] = dict()
        cls._logs[EXEC_MSG] = dict()
        cls._logs[EXEC_ERROR] = dict()
        cls._logs[DB_MERGE_RESULT] = list()
        cls._dbs_in_log = set()
        cls._tagged_schema_db_name = None

    '''Initializes all logs for a single database (private function) if they do not exist'''
    @classmethod
    def _init_logs_for_new_db(cls, db_name):
        if db_name not in cls._dbs_in_log:
            cls._logs[INVALID_ACC_NUM_FMT][db_name] = list()
            cls._logs[ACC_NUM_NOT_FOUND][db_name] = list()
            cls._logs[FALSE_SEQ_CLASS][db_name] = list()
            cls._logs[UNKNOWN_SEQ_TYPE][db_name] = list()
            cls._logs[SEQ_MISMATCH][db_name] = list()
            cls._logs[DUPLICATED_HEADER][db_name] = list()
            cls._logs[REDUNDANT_SEQ][db_name] = list()
            cls._logs[OBSOLETE_VER][db_name] = list()
            cls._logs[EXEC_MSG][db_name] = list()
            cls._logs[EXEC_ERROR][db_name] = list()
            cls._dbs_in_log.add(db_name)

    '''
    Function name: _tag_schema_db
    Inputs       : Database name
    Outputs      : Database name with/without tag
    Description  : Tags the input database name with the schema database label SCHEMA_DB_TAG if the
                   input database is not the default one
    '''
    @classmethod
    def _tag_schema_db(cls, db_name):
        if db_name != DEFAULT_DB_NAME:
            cls._tagged_schema_db_name = '{}{}'.format(db_name, SCHEMA_DB_TAG)
            return cls._tagged_schema_db_name

        return db_name

    '''
    Function name: _log_msg
    Inputs       : Message, log name, database name, boolean indicating schema database
    Outputs      : Nil
    Description  : Base function to write input message to a specific log for a specific database
    '''
    @classmethod
    def _log_msg(cls, msg, log_name, db_name = DEFAULT_DB_NAME, is_for_schema_db = False):
        if is_for_schema_db:
            db_name = cls._tag_schema_db(db_name)

        cls._init_logs_for_new_db(db_name)
        cls._logs[log_name][db_name].append('{}{}'.format(msg, os.linesep))

    '''Log for invalid NCBI accession number format'''
    log_invalid_acc_num_fmt = partialmethod(_log_msg, log_name = INVALID_ACC_NUM_FMT)
    '''Log for non-existing NCBI accession number'''
    log_acc_num_not_found = partialmethod(_log_msg, log_name = ACC_NUM_NOT_FOUND)
    '''log_cds_not_found = partialmethod(_log_msg, log_name = CDS_NOT_FOUND)'''
    '''Log for false sequence class annotation'''
    log_false_seq_class = partialmethod(_log_msg, log_name = FALSE_SEQ_CLASS)
    '''Log for unknown input sequence type'''
    log_unknown_seq_type = partialmethod(_log_msg, log_name = UNKNOWN_SEQ_TYPE)
    '''log_mixed_seq_type = partialmethod(_log_msg, log_name = MIXED_SEQ_TYPE)'''
    '''Log for duplicated sequence header for two sequences'''
    log_duplicated_header = partialmethod(_log_msg, log_name = DUPLICATED_HEADER)
    '''Log for obsolete NCBI accession version for sequence'''
    log_obsolete_ver = partialmethod(_log_msg, log_name = OBSOLETE_VER)
    '''Log for execution message'''
    log_exec_msg = partialmethod(_log_msg, log_name = EXEC_MSG)
    '''Log for execution error'''
    log_exec_error = partialmethod(_log_msg, log_name = EXEC_ERROR)
    '''Log for data retrieval error'''
    log_data_retrieval_error = partialmethod(_log_msg, msg = DATA_RETRIEVAL_ERROR, log_name = EXEC_ERROR)

    '''
    Function name: log_seq_mismatch
    Inputs       : Sequence Id, obsolete NCBI accession version flag, database name, boolean indicating
                   schema database
    Outputs      : Nil
    Description  : Log for sequence mismatch for a specific database
    '''
    @classmethod
    def log_seq_mismatch(cls, seq_id, is_obsolete_ver, db_name = DEFAULT_DB_NAME, is_for_schema_db = False):
        if is_for_schema_db:
            db_name = cls._tag_schema_db(db_name)

        cls._init_logs_for_new_db(db_name)
        if is_obsolete_ver:
            cls._logs[SEQ_MISMATCH][db_name].append('{} (Note: Obsolete sequence){}'.format(seq_id, os.linesep))
        else:
            cls._logs[SEQ_MISMATCH][db_name].append('{}{}'.format(seq_id, os.linesep))

    '''
    Function name: log_redundant_seq
    Inputs       : Sequence Id, sequence Id with redundant sequence, reverse complementarity flag,
                   database name, boolean indicating schema database
    Outputs      : Nil
    Description  : Log for redundant sequence for a specific database
    '''
    @classmethod
    def log_redundant_seq(cls, seq_id, redundant_seq_id, is_rev_comp, db_name = DEFAULT_DB_NAME,
                          is_for_schema_db = False):
        if is_for_schema_db:
            db_name = cls._tag_schema_db(db_name)

        cls._init_logs_for_new_db(db_name)
        if is_rev_comp:
            redundant_seq_msg = '{} (reverse complement of {}){}'.format(seq_id, redundant_seq_id, os.linesep)
        else:
            redundant_seq_msg = '{} (redundant with {}){}'.format(seq_id, redundant_seq_id, os.linesep)

        cls._logs[REDUNDANT_SEQ][db_name].append(redundant_seq_msg)

    '''
    Function name: _export_log
    Inputs       : Output stream, boolean controlling the export of empty log, log name to export,
                   header of the log
    Outputs      : Nil
    Description  : Base function to export a specific log for all databases to the output stream
    '''
    @classmethod
    def _export_log(cls, output_stream, is_export_empty_log, log_name, log_header):
        msg_log = list()

        if cls._tagged_schema_db_name is not None:
            db_msg_log = cls._logs[log_name][cls._tagged_schema_db_name]
            if len(db_msg_log) > 0:
                db_name = cls._tagged_schema_db_name.replace(SCHEMA_DB_TAG, '')
                msg_log.append('For schema database {}:{}'.format(db_name, os.linesep))
                db_msg_log.sort()
                msg_log += db_msg_log
                msg_log.append(os.linesep)
        
        is_print_db_name = (len(cls._dbs_in_log) > 1)

        for db_name in sorted(cls._dbs_in_log):
            db_msg_log = cls._logs[log_name][db_name]
            if len(db_msg_log) == 0 or db_name == cls._tagged_schema_db_name:
                continue

            if is_print_db_name:
                msg_log.append('For {}:{}'.format(db_name, os.linesep))
                
            db_msg_log.sort()
            msg_log += db_msg_log
            msg_log.append(os.linesep)

        if len(msg_log) == 0:
            if is_export_empty_log:
                output_stream.writelines('{}{}Nil{}{}'.format(log_header, os.linesep, os.linesep, os.linesep))
        else:
            output_stream.writelines('{}{}'.format(log_header, os.linesep))
            output_stream.writelines(msg_log)
            '''output_stream.writelines(os.linesep)'''

    '''Export log for invalid NCBI accession number format'''
    _export_invalid_acc_num_fmt_log = partialmethod(_export_log, log_name = INVALID_ACC_NUM_FMT,
                                                    log_header = '----- No valid accession number format -----')
    '''Export log for non-existing NCBI accession number'''
    _export_acc_num_not_found_log = partialmethod(_export_log, log_name = ACC_NUM_NOT_FOUND,
                                                  log_header = '----- Accession numbers removed/not found -----')
    '''
    _export_cds_not_found_log = partialmethod(_export_log, log_name = CDS_NOT_FOUND,
                                              log_header = '----- Matching CDS not found -----')
    '''
    '''Export log for false sequence class annotation'''
    _export_false_seq_class_log = partialmethod(_export_log, log_name = FALSE_SEQ_CLASS,
                                                log_header = '----- Potential incorrect sequence classification -----')
    '''Export log for unknown input sequence type'''
    _export_unknown_seq_type_log = partialmethod(_export_log, log_name = UNKNOWN_SEQ_TYPE,
                                                 log_header = '----- Unknown sequence type -----')
    '''
    _export_mixed_seq_type_log = partialmethod(_export_log, log_name = MIXED_SEQ_TYPE,
                                               log_header = '----- Inconsistent base sequence type -----')
    '''
    '''Export log for sequence mismatch'''
    _export_seq_mismatch_log = partialmethod(_export_log, log_name = SEQ_MISMATCH,
                                             log_header = '----- Sequence mis-matches -----')
    '''Export log for duplicated sequence header'''
    _export_duplicated_header_log = partialmethod(_export_log, log_name = DUPLICATED_HEADER,
                                                  log_header = '----- Duplicated sequence headers -----')
    '''Export log for redundant sequence'''
    _export_redundant_seq_log = partialmethod(_export_log, log_name = REDUNDANT_SEQ,
                                              log_header = '----- Redundant sequences -----')
    '''Export log for obsolete NCBI accession version for sequence'''
    _export_obsolete_ver_log = partialmethod(_export_log, log_name = OBSOLETE_VER,
                                             log_header = '----- Obsolete sequence version -----')
    '''Export log for execution message'''
    export_exec_msg = partialmethod(_export_log, log_name = EXEC_MSG, log_header = '----- Messages -----')
    '''Export log for execution error'''
    export_exec_error = partialmethod(_export_log, is_export_empty_log = False, log_name = EXEC_ERROR,
                                      log_header = '----- Execution errors -----')

    '''
    Function name: export_qc_check_logs
    Inputs       : Output stream, boolean controlling the check of sequence classfication
    Outputs      : Nil
    Description  : Export all logs for ARG database validation, including empty log, to the output
                   stream
    '''
    @classmethod
    def export_qc_check_logs(cls, output_stream, is_check_seq_class):
        cls._export_invalid_acc_num_fmt_log(output_stream, is_export_empty_log = True)
        cls._export_acc_num_not_found_log(output_stream, is_export_empty_log = True)
        '''cls._export_cds_not_found_log(output_stream, has_nt_seq_check)'''
        cls._export_false_seq_class_log(output_stream, is_export_empty_log = is_check_seq_class)
        cls._export_unknown_seq_type_log(output_stream, is_export_empty_log = True)
        '''cls._export_mixed_seq_type_log(output_stream, True)'''
        cls._export_seq_mismatch_log(output_stream, is_export_empty_log = True)
        cls._export_duplicated_header_log(output_stream, is_export_empty_log = True)
        cls._export_redundant_seq_log(output_stream, is_export_empty_log = True)
        cls._export_obsolete_ver_log(output_stream, is_export_empty_log = True)
        cls.export_exec_msg(output_stream, is_export_empty_log = False)

    '''
    Function name: export_qc_check_summary
    Inputs       : Output stream, total number of sequences, boolean controlling the check of
                   sequence class annotation, external summary, database name
    Outputs      : Nil
    Description  : Export summary for all logs for a specific database to the output stream
    '''
    @classmethod
    def export_qc_check_summary(cls, output_stream, total_seq_record_count, is_check_seq_class,
                                ext_summary = None, db_name = DEFAULT_DB_NAME):
        cls._init_logs_for_new_db(db_name)
        summary = list()
        summary_stmt = cls.create_summary_stmt(total_seq_record_count, 'inspected:', num_display_width = 0)
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(len(cls._logs[INVALID_ACC_NUM_FMT][db_name]),
                                               'with no valid accession number format')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(len(cls._logs[ACC_NUM_NOT_FOUND][db_name]),
                                               'with accession number removed/not found')
        summary.append(summary_stmt)

        '''
        if has_nt_seq_check:
            summary_stmt = cls.create_summary_stmt(len(cls._logs[CDS_NOT_FOUND]), 'with no matching CDS found')
            summary.append(summary_stmt)
        '''

        if is_check_seq_class:
            summary_stmt = cls.create_summary_stmt(len(cls._logs[FALSE_SEQ_CLASS][db_name]),
                                                   'with potential incorrect sequence classification')
            summary.append(summary_stmt)

        summary_stmt = cls.create_summary_stmt(len(cls._logs[UNKNOWN_SEQ_TYPE][db_name]), 'with unknown sequence type')
        summary.append(summary_stmt)
        '''
        summary_stmt = cls.create_summary_stmt(len(cls._logs[MIXED_SEQ_TYPE]), 'with mixed base sequence type')
        summary.append(summary_stmt)
        '''
        summary_stmt = cls.create_summary_stmt(len(cls._logs[SEQ_MISMATCH][db_name]), 'with sequence mis-match')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(len(cls._logs[DUPLICATED_HEADER][db_name]), 'with duplicated headers')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(len(cls._logs[REDUNDANT_SEQ][db_name]), 'being redundant')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(len(cls._logs[OBSOLETE_VER][db_name]), 'obsolete')
        summary.append(summary_stmt)

        output_stream.writelines('----- Summary -----{}'.format(os.linesep))
        if ext_summary is None:
            output_stream.writelines(summary)
        else:
            output_stream.writelines(summary + ext_summary)

    '''
    Function name: export_ext_summary
    Inputs       : Output stream, external summary
    Outputs      : Nil
    Description  : Export external summary to the output stream
    '''
    @staticmethod
    def export_ext_summary(output_stream, ext_summary):
        output_stream.writelines('----- Summary -----{}'.format(os.linesep))
        output_stream.writelines(ext_summary)

    '''
    Function name: create_summary_stmt
    Inputs       : Number of sequences, suffix of the summary statement, description of the sequence,
                   width for the numeric display field
    Outputs      : Summary statement
    Description  : Generate summary statement
    '''
    @staticmethod
    def create_summary_stmt(seq_record_count, stmt_suffix, seq_desc = None, num_display_width = 5):
        if seq_record_count < 2:
            seq_term = 'sequence'
        else:
            seq_term = 'sequences'

        if seq_desc is None:
            return '{} {} {}{}'.format(str.rjust(str(seq_record_count), num_display_width), seq_term, stmt_suffix,
                                       os.linesep)
        else:
            return '{} {} {} {}{}'.format(str.rjust(str(seq_record_count), num_display_width), seq_desc, seq_term,
                                          stmt_suffix, os.linesep)

    '''
    Function name: has_exec_error
    Inputs       : Database name
    Outputs      : Boolean indicating occurrence of execution error 
    Description  : Check whether execution error has occurred for a particular database
    '''
    @classmethod
    def has_exec_error(cls, db_name = DEFAULT_DB_NAME):
        cls._init_logs_for_new_db(db_name)
        for db_error_msg_log in cls._logs[EXEC_ERROR].values():
            if len(db_error_msg_log) > 0:
                return True

        return False

    '''
    Function name: log_merge_db_result
    Inputs       : Header of the exported sequence, and the mapping between this exported
                   sequence and databases containing it
    Outputs      : Nil
    Description  : Log for source databases contributing the exported sequence
    '''
    @classmethod
    def log_merge_db_result(cls, exported_seq_header, redundant_seq_rec_src_tags):
        merge_output_seq_msg = cls.MERGE_OUTPUT_SEQ_SRC_MSG_TEMPLATE.format(exported_seq_header, os.linesep)
        cls._logs[DB_MERGE_RESULT].append(merge_output_seq_msg)

        for seq_record_src_tag in redundant_seq_rec_src_tags:
            src_db_name = os.path.basename(seq_record_src_tag[1])
            seq_src_db_msg = cls.SEQ_SRC_DB_MSG_TEMPLATE.format(src_db_name, seq_record_src_tag[2], os.linesep)
            cls._logs[DB_MERGE_RESULT].append(seq_src_db_msg)

        if hasattr(cls, '_redundant_seq_count'):
            cls._redundant_seq_count += len(redundant_seq_rec_src_tags) - 1
        else:
            cls._redundant_seq_count = len(redundant_seq_rec_src_tags) - 1

    '''
    Function name: log_cross_db_duplicated_header
    Inputs       : Mapping between sequence with duplicated (auto-generated) header, database
                   containing this sequence, and the original sequence header (if applicable)
    Outputs      : Nil
    Description  : Log for duplicated sequence header across two databases
    '''
    @classmethod
    def log_cross_db_duplicated_header(cls, seq_record_src_tag, is_auto_annotate):
        '''If the auto-generated annotation header is duplicated, show also the original header'''
        if is_auto_annotate:
            autogen_seq_header = seq_record_src_tag[0].id.replace(SEQ_VERSION_MARKUP, '')
            exec_msg = cls.CROSS_DB_DUPLICATED_AUTOGEN_HDR_MSG_TEMPLATE.format(autogen_seq_header,
                                                                               seq_record_src_tag[2])
        else:
            exec_msg = cls.CROSS_DB_DUPLICATED_SEQ_HDR_MSG_TEMPLATE.format(seq_record_src_tag[2])

        cls.log_exec_msg(exec_msg, db_name = seq_record_src_tag[1])
        if hasattr(cls, '_duplicated_hdr_count'):
            cls._duplicated_hdr_count += 1
        else:
            cls._duplicated_hdr_count = 1

    '''
    Function name: _export_merge_db_result
    Inputs       : Output stream
    Outputs      : Nil
    Description  : Export log for the sequences generated in database merge
    '''
    @classmethod
    def _export_merge_db_result(cls, output_stream):
        output_stream.writelines('{}{}'.format('----- Database consolidation results -----', os.linesep))

        if len(cls._logs[DB_MERGE_RESULT]) == 0:
            output_stream.writelines('{}{}{}'.format('Nil', os.linesep, os.linesep))
        else:
            output_stream.writelines(cls._logs[DB_MERGE_RESULT])

    '''
    Function name: export_merge_db_check_logs
    Inputs       : Output stream
    Outputs      : Nil
    Description  : Export all relevant logs for ARG database integration to the output stream
    '''
    @classmethod
    def export_merge_db_check_logs(cls, output_stream):
        cls._export_invalid_acc_num_fmt_log(output_stream, is_export_empty_log = False)
        cls._export_acc_num_not_found_log(output_stream, is_export_empty_log = False)
        cls._export_unknown_seq_type_log(output_stream, is_export_empty_log = False)
        cls._export_seq_mismatch_log(output_stream, is_export_empty_log = False)
        cls.export_exec_msg(output_stream, is_export_empty_log = False)
        cls._export_merge_db_result(output_stream)

    '''
    Function name: export_merge_db_summary
    Inputs       : Output stream, number of source database sequences, number of schema database
                   sequences, external summary
    Outputs      : Nil
    Description  : Export summary for all relevant logs for all databases to the output stream
    '''
    @classmethod
    def export_merge_db_summary(cls, output_stream, src_db_seq_rec_count, schema_db_seq_rec_count = 0,
                                ext_summary = None):
        summary = list()

        if schema_db_seq_rec_count > 0 and cls._tagged_schema_db_name is not None:
            db_name = cls._tagged_schema_db_name.replace(SCHEMA_DB_TAG, '')
            summary_stmt_suffix = 'read from schema database {}:'.format(db_name)
            summary_stmt = cls.create_summary_stmt(schema_db_seq_rec_count, summary_stmt_suffix, num_display_width = 0)
            summary.append(summary_stmt)
            summary_stmt = cls.create_summary_stmt(len(cls._logs[INVALID_ACC_NUM_FMT][cls._tagged_schema_db_name]),
                                                   'with no valid accession number format, skipped')
            summary.append(summary_stmt)
            summary_stmt = cls.create_summary_stmt(len(cls._logs[ACC_NUM_NOT_FOUND][cls._tagged_schema_db_name]),
                                                   'with accession number removed/not found, skipped')
            summary.append(summary_stmt)
            summary_stmt = cls.create_summary_stmt(len(cls._logs[UNKNOWN_SEQ_TYPE][cls._tagged_schema_db_name]),
                                                   'with unknown sequence type, skipped')
            summary.append(summary_stmt)
            summary_stmt = cls.create_summary_stmt(len(cls._logs[SEQ_MISMATCH][cls._tagged_schema_db_name]),
                                                   'with sequence mis-match, skipped')
            summary.append(summary_stmt)
            summary.append(os.linesep)

        src_db_invalid_acc_num_fmt_count = 0
        src_db_acc_num_not_found_count = 0
        src_db_unknown_seq_type_count = 0
        src_db_seq_mismatch_count = 0

        for db_name in cls._dbs_in_log:
            if cls._tagged_schema_db_name is None or db_name != cls._tagged_schema_db_name:
                src_db_invalid_acc_num_fmt_count += len(cls._logs[INVALID_ACC_NUM_FMT][db_name])
                src_db_acc_num_not_found_count += len(cls._logs[ACC_NUM_NOT_FOUND][db_name])
                src_db_unknown_seq_type_count += len(cls._logs[UNKNOWN_SEQ_TYPE][db_name])
                src_db_seq_mismatch_count += len(cls._logs[SEQ_MISMATCH][db_name])

        summary_stmt = cls.create_summary_stmt(src_db_seq_rec_count, 'read from source databases:',
                                               num_display_width = 0)
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(src_db_invalid_acc_num_fmt_count,
                                               'with no valid accession number format, skipped')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(src_db_acc_num_not_found_count,
                                               'with accession number removed/not found, skipped')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(src_db_unknown_seq_type_count, 'with unknown sequence type, skipped')
        summary.append(summary_stmt)
        summary_stmt = cls.create_summary_stmt(src_db_seq_mismatch_count, 'with sequence mis-match, skipped')
        summary.append(summary_stmt)

        if hasattr(cls, '_redundant_seq_count'):
            redundant_seq_stmt = ProcLog.create_summary_stmt(cls._redundant_seq_count, 'found, skipped', 'redundant')
        else:
            redundant_seq_stmt = ProcLog.create_summary_stmt(0, 'found, skipped', 'redundant')
            
        summary.append(redundant_seq_stmt)

        if hasattr(cls, '_duplicated_hdr_count'):
            duplicated_hdr_stmt = ProcLog.create_summary_stmt(cls._duplicated_hdr_count,
                                                              'with duplicated header, skipped')
        else:
            duplicated_hdr_stmt = ProcLog.create_summary_stmt(0, 'with duplicated header, skipped')

        summary.append(duplicated_hdr_stmt)

        output_stream.writelines('----- Summary -----{}'.format(os.linesep))
        if ext_summary is None:
            output_stream.writelines(summary)
        else:
            output_stream.writelines(summary + ext_summary)
