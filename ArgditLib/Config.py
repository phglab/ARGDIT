from .ProcLog import ProcLog
import configparser
import re

ARDGIT_SECTION = 'ARGDIT'
FASTA_HEADER_FIELD_SEP = 'FastaHeaderFieldSeparator'
FIELD_SEP = 'OperationalFieldSeparator'
OTL_ANNOT_CHECK_SECTION = 'Ontology annotation check'
MIN_SEQ_COUNT = 'MinSequenceCount'
BOOTSTRAP_FACTOR = 'BootstrapFactor'
ENTREZ_SECTION = 'Entrez'
EMAIL = 'Email'
EMAIL_PATTERN = r'(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)'

'''Object class for configuration file'''
class Config:
    def __init__(self, config_file_path):
        config_settings = configparser.ConfigParser()
        config_settings.read(config_file_path)

        try:
            self._fasta_header_field_sep = config_settings[ARDGIT_SECTION][FASTA_HEADER_FIELD_SEP]
            self._field_sep = config_settings[ARDGIT_SECTION][FIELD_SEP]
            self._min_seq_count = int(config_settings[OTL_ANNOT_CHECK_SECTION][MIN_SEQ_COUNT])
            self._bootstrap_factor = int(config_settings[OTL_ANNOT_CHECK_SECTION][BOOTSTRAP_FACTOR])
            self._entrez_email = config_settings[ENTREZ_SECTION][EMAIL]
        except KeyError as e:
            ProcLog.log_exec_error('Configuration not set: {}'.format(str(e)))

        if not re.match(EMAIL_PATTERN, self._entrez_email):
            ProcLog.log_exec_error('Invalid email address: {}'.format(self._entrez_email))

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

