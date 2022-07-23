'''Constants referred by the ARGDIT tools'''

REFSEQ_NT_PREFIX_PATTERN = r'NC_|AC_|NG_|NT_|NW_|NZ_|NM_|XM_'
WGS_ACC_NUM_FMT_LIST = [r'((' + REFSEQ_NT_PREFIX_PATTERN + ')?[A-Z]{4}\d{8,}(\.\d)?)',
                        r'((' + REFSEQ_NT_PREFIX_PATTERN + ')?[A-Z]{6}\d{9,}(\.\d)?)']
WGS_ACC_NUM_PATTERN = r'({})'.format(r'|'.join(WGS_ACC_NUM_FMT_LIST))
NT_ACC_NUM_FMT_LIST = [r'((' + REFSEQ_NT_PREFIX_PATTERN + ')?[A-Z]\d{5}(\.\d)?)',
                       r'((' + REFSEQ_NT_PREFIX_PATTERN + ')?[A-Z]{2}\d{6}(\.\d)?)',
                       r'((' + REFSEQ_NT_PREFIX_PATTERN + ')?[A-Z]{2}\d{8}(\.\d)?)',
                       r'((' + REFSEQ_NT_PREFIX_PATTERN + ')\d{6,}(\.\d)?)']
NT_ACC_NUM_PATTERN = r'({}|{})'.format(r'|'.join(WGS_ACC_NUM_FMT_LIST), r'|'.join(NT_ACC_NUM_FMT_LIST))

PROTEIN_ACC_NUM_FMT_LIST = [r'([A-Z]{3}\d{5}(\.\d)?)',
                            r'([A-Z]{3}\d{7}(\.\d)?)',
                            r'((WP_|AP_|NP_|YP_|XP_)\d+(\.\d)?)']
PROTEIN_ACC_NUM_PATTERN = r'({})'.format(r'|'.join(PROTEIN_ACC_NUM_FMT_LIST))

NT_ACC_NUM_EMBED_PATTERN = r'([^A-Z]|^){}'.format(NT_ACC_NUM_PATTERN)
PROTEIN_ACC_NUM_EMBED_PATTERN = r'([^A-Z]|^){}'.format(PROTEIN_ACC_NUM_PATTERN)

IUPAC_DNA_STR_PATTERN = r'^[ACGTRYSWKMBDHVN]+$'
IUPAC_PROTEIN_STR_PATTERN = r'^[ABCDEFGHIKLMNPQRSTVWXYZ]+$'
IUPAC_AMBIG_DNA_BASES = r'[RYSWKMBDHVN]'
IUPAC_AMBIG_PROTEIN_BASES = r'[X]'

SEQ_VERSION_MARKUP = '<ARGDIT seq version>'
NCBI_LIVE_SEQ_STATUS = 'live'
