from .ProcLog import ProcLog
import configparser
import re

ARDGIT_SECTION = 'ARGDIT'
FASTA_HEADER_FIELD_SEP = 'FastaHeaderFieldSeparator'
FIELD_SEP = 'OperationalFieldSeparator'
SEQ_CLASS_CHECK_SECTION = 'Sequence classification check'
MIN_SEQ_COUNT = 'MinSequenceCount'
BOOTSTRAP_FACTOR = 'BootstrapFactor'
ENTREZ_SECTION = 'Entrez'
EMAIL = 'Email'
TRANSLATION_SECTION = 'Translate'
DEFAULT_GENETIC_CODE = 'DefaultGeneticCode'
EMAIL_PATTERN = r'(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)'

'''Object class for configuration file'''
class Config:
    def __init__(self, config_file_path):
        config_settings = configparser.ConfigParser()
        config_settings.read(config_file_path)

        try:
            self._fasta_header_field_sep = config_settings[ARDGIT_SECTION][FASTA_HEADER_FIELD_SEP]
            self._field_sep = config_settings[ARDGIT_SECTION][FIELD_SEP]
            self._min_seq_count = config_settings[SEQ_CLASS_CHECK_SECTION][MIN_SEQ_COUNT]
            self._bootstrap_factor = config_settings[SEQ_CLASS_CHECK_SECTION][BOOTSTRAP_FACTOR]
            self._entrez_email = config_settings[ENTREZ_SECTION][EMAIL]
            self._default_genetic_code = config_settings[TRANSLATION_SECTION][DEFAULT_GENETIC_CODE]
        except KeyError as e:
            ProcLog.log_exec_error('Configuration not set: {}'.format(str(e)))

        if not re.match(EMAIL_PATTERN, self._entrez_email):
            ProcLog.log_exec_error('Invalid email address format {} in configuration file'.format(self._entrez_email))

        if not re.match('\d+', self._min_seq_count):
            err_msg = 'Non-integer minimum sequence count {} in configuration file'.format(self._min_seq_count)
            ProcLog.log_exec_error(err_msg)

        if not re.match('\d+', self._bootstrap_factor):
            err_msg = 'Non-integer bootstrap factor {} in configuration file'.format(self._bootstrap_factor)
            ProcLog.log_exec_error(err_msg)

        if not re.match('\d+', self._default_genetic_code):
            err_msg = 'Non-integer default genetic code {} in configuration file'.format(self._default_genetic_code)
            ProcLog.log_exec_error(err_msg)

        self._min_seq_count = int(self._min_seq_count)
        self._bootstrap_factor = int(self._bootstrap_factor)
        self._default_genetic_code = int(self._default_genetic_code)

    @property
    def fasta_header_field_sep(self):
        return self._fasta_header_field_sep

    @property
    def field_sep(self):
        return self._field_sep

    @property
    def min_seq_count(self):
        return self._min_seq_count

    @property
    def bootstrap_factor(self):
        return self._bootstrap_factor

    @property
    def entrez_email(self):
        return self._entrez_email

    @property
    def default_genetic_code(self):
        return self._default_genetic_code
