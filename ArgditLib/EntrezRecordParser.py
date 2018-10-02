'''Collection of data handle parser classes NCBI database searching'''

from .CDSRegion import CDSRegion
from .Constants import NT_ACC_NUM_PATTERN, PROTEIN_ACC_NUM_PATTERN
from .ProteinInfo import ProteinInfo
from .Utils import extract_protein_acc_num, trim_version
from Bio import Entrez
from Bio import SeqIO
from xml.etree import ElementTree
import re

HEADER_NT_ACC_NUM_PATTERN = r'^>Feature [a-z]+\|' + NT_ACC_NUM_PATTERN + r'\|.*\n$'
HEADER_PROTEIN_ACC_NUM_PATTERN = r'^>Feature [a-z]+\|' + PROTEIN_ACC_NUM_PATTERN + r'\|.*\n$'
CDS_PROTEIN_ACC_NUM_PATTERN = r'^\t+protein_id.*\|' + PROTEIN_ACC_NUM_PATTERN + r'\|.*\n$'
GENE_SHORT_NAME_PATTERN = r'^\t+gene\t+(.+?)\n$'
PRODUCT_PATTERN = r'^\t+product\t+(.+?)\n$'
TRANSLATE_EXCEPT_PATTERN = r'^\t+transl_except\t+\(pos:(\d+)\.\.(\d+),aa:Met\)\n$'
CODON_START_PATTERN = r'^\t+codon_start\t+(\d+)\n$'
ANY_REGION_PATTERN = r'^(<?)(\d+)\t(>?)(\d+)\t.+\n$'
GENE_REGION_PATTERN = r'^(<?)(\d+)\t(>?)(\d+)\tgene\n$'
CDS_REGION_PATTERN = r'^(<?)(\d+)\t(>?)(\d+)\tCDS\n$'
PROTEIN_REGION_PATTERN = r'^(<?)(\d+)\.\.(>?)(\d+)$'
ADDITIONAL_REGION_PATTERN = r'^(<?)(\d+)\t(>?)(\d+)\n$'
ERROR_PATTERN = r'^\t<ERROR>.+</ERROR>\n$'

'''EPost parser class to extract the query key and web environment returned from Entrez epost function'''
class EPostParser:
    '''
    Function name: parse
    Inputs       : Data handle
    Outputs      : Query key, web environment, and errors returned, if any
    Description  : Parses the data handle returned from the Entrez epost utility. The data is in XML
                   format
    '''
    def parse(handle):
        query_key = None
        web_env = None
        errors = None
        
        xml_tree = ElementTree.parse(handle)
        root_node = xml_tree.getroot()
        
        for child_node in root_node:
            if child_node.tag == 'QueryKey':
                query_key = child_node.text
            elif child_node.tag == 'WebEnv':
                web_env = child_node.text
            elif child_node.tag == 'ERROR':
                if errors is None:
                    errors = [child_node.text]
                else:
                    errors.append(child_node.text)

        return query_key, web_env, errors

'''
Feature table parser class to extract coding region (CDS) from the data handle returned from the
Entrez efetch function

---Attributes---
cds_seq_len_filters: CDS sequence length filters
target_cds_region_grps: Target CDS information matching respective sequence length filter, if any.
                        Lists of CDSRegion objects categorized according to their associated nucleotide
                        accession numbers
target_protein_ids: Accession numbers (or IDs) of the target proteins
is_parse_complete: Boolean indicating whether the parse is completed
'''
class FeatTblCDSParser:
    '''Constructor'''
    def __init__(self, cds_seq_len_filters = None):
        self._cds_seq_len_filters = cds_seq_len_filters
        self._target_cds_region_grps = dict()
        self._target_protein_ids = set()
        self._is_parse_complete = True

    '''
    def get_protein_to_nt_acc_num_map(self):
        protein_to_nt_acc_num_map = dict()

        for nt_acc_num, target_cds_region_grps in self._target_cds_region_grps.items():
            for target_cds in target_cds_region_grps:
                if target_cds.protein_id in protein_to_nt_acc_num_map:
                    protein_to_nt_acc_num_map[target_cds.protein_id].append(nt_acc_num)
                else:
                    protein_to_nt_acc_num_map[target_cds.protein_id] = [nt_acc_num]

        return protein_to_nt_acc_num_map
    '''

    def get_target_cds_region_groups(self):
        return self._target_cds_region_grps

    def get_target_protein_ids(self):
        return self._target_protein_ids

    @property
    def is_parse_complete(self):
        return self._is_parse_complete

    '''
    Function name: _is_target_cds
    Inputs       : CDS information, CDS sequence length filter
    Outputs      : Boolean
    Description  : Determines whether the CDS region matches the target CDS sequence length
    '''
    def _is_target_cds(self, cds_region, cds_seq_len_filter):
        if cds_seq_len_filter is None:
            return True

        return cds_region.length in cds_seq_len_filter

    '''
    Function name: _update_target_cds
    Inputs       : Nucleotide accession number, CDS information, and CDS sequence length filter
    Outputs      : Nil
    Description  : When the input CDS region matches the target CDS sequence length, adds it to the
                   target CDS information and records its protein accession number (ID)
    '''
    def _update_target_cds(self, nt_acc_num, cds_region, cds_seq_len_filter):
        if self._is_target_cds(cds_region, cds_seq_len_filter):
            self._target_cds_region_grps[nt_acc_num].append(cds_region)
            self._target_protein_ids.add(cds_region.protein_id)

    '''
    Function name: parse
    Inputs       : Feature table data handle
    Outputs      : None
    Description  : Parses the data handle to extract all CDS information, and retains those matching
                   the CDS sequence length filter, if any
    '''
    def parse(self, ft_handle):       
        nt_acc_num = None
        cds_region = CDSRegion()
        cds_seq_len_filter = None
        gene_region_range = None
        gene_short_name = None
        is_in_cds_region = False
        is_last_line_define_cds_region = False

        while True:
            line = ft_handle.readline()
            '''An empty line indicates the end of data'''
            if line == '':
                if nt_acc_num is not None and cds_region.protein_id is not None:
                    self._update_target_cds(nt_acc_num, cds_region, cds_seq_len_filter)
                break

            if re.match(ERROR_PATTERN, line):
                self._is_parse_complete = False
                break

            '''Beginning of a new nucleotide accession number'''
            m = re.match(HEADER_NT_ACC_NUM_PATTERN, line)
            if m:
                '''Save the last parsed CDS information'''
                if nt_acc_num is not None and cds_region.protein_id is not None:
                    self._update_target_cds(nt_acc_num, cds_region, cds_seq_len_filter)

                nt_acc_num = m.group(1)
                cds_region = CDSRegion()
                self._target_cds_region_grps[nt_acc_num] = []
                '''Load the corresponding CDS sequence length filter'''
                if self._cds_seq_len_filters is not None:
                    nt_non_ver_acc_num = trim_version(nt_acc_num)
                    if nt_non_ver_acc_num in self._cds_seq_len_filters:
                        cds_seq_len_filter = self._cds_seq_len_filters[nt_non_ver_acc_num]
                    else:
                        cds_seq_len_filter = set()

                gene_region_range = None
                is_in_cds_region = False
                continue

            if re.match(ANY_REGION_PATTERN, line):
                '''Beginning of a new CDS region'''
                m = re.match(CDS_REGION_PATTERN, line)
                if m:
                    '''Save the last parsed CDS information'''
                    if nt_acc_num is not None and cds_region.protein_id is not None:
                        self._update_target_cds(nt_acc_num, cds_region, cds_seq_len_filter)

                    cds_region = CDSRegion()
                    '''Sequence completeness attribute in CDS information'''
                    cds_region.add_region_range(int(m.group(2)), int(m.group(4)))
                    if m.group(1) == '<':
                        cds_region.is_5_partial = True
                    if m.group(3) == '>':
                        cds_region.is_3_partial = True

                    is_in_cds_region = True
                    is_last_line_define_cds_region = True
                    continue
                else:
                    is_in_cds_region = False

                '''
                Beginning of a new gene region which precedes its associated CDS region, if this CDS
                region exists
                '''
                m = re.match(GENE_REGION_PATTERN, line)
                if m:
                    gene_region_range = (int(m.group(2)), int(m.group(4)))
                    gene_short_name = None

                continue

            '''
            If a single CDS region consists of multiple region ranges, matches the rest ranges after
            the first range matched above
            '''
            m = re.match(ADDITIONAL_REGION_PATTERN, line)
            if m and is_in_cds_region:
                cds_region.add_region_range(int(m.group(2)), int(m.group(4)))
                if m.group(1) == '<':
                    cds_region.is_5_partial = True
                if m.group(3) == '>':
                    cds_region.is_3_partial = True

                continue

            '''
            Upon the CDS region range definition finishes (and so the entire CDS range is known), when
            the its gene region is also defined, put also the abbreviated gene name in the CDS
            information
            '''
            if is_last_line_define_cds_region:
                if gene_region_range is not None and gene_short_name is not None:
                    cds_region_range = cds_region.get_region_range()
                    if cds_region.is_complementary:
                        if cds_region_range[0] <= gene_region_range[0] and cds_region_range[1] >= gene_region_range[1]:
                            cds_region.gene_short_name = gene_short_name
                    else:
                        if cds_region_range[0] >= gene_region_range[0] and cds_region_range[1] <= gene_region_range[1]:
                            cds_region.gene_short_name = gene_short_name

                is_last_line_define_cds_region = False

            '''Codon start attribute in CDS information'''
            m = re.match(CODON_START_PATTERN, line)
            if m and is_in_cds_region:
                cds_region.codon_start = int(m.group(1))
                continue

            '''Special start codon attribute in CDS information'''
            m = re.match(TRANSLATE_EXCEPT_PATTERN, line)
            if m and is_in_cds_region:
                except_codon_start = int(m.group(1))
                if abs(int(m.group(2)) - except_codon_start) == 2:
                    cds_region.set_enforce_start_codon(except_codon_start)

                continue

            '''Protein accession number attribute in CDS information'''
            m = re.match(CDS_PROTEIN_ACC_NUM_PATTERN, line)
            if m and is_in_cds_region:
                cds_region.protein_id = m.group(1)
                continue

            '''Abbreviated gene name attribute in CDS information'''
            m = re.match(GENE_SHORT_NAME_PATTERN, line)
            if m and gene_region_range is not None:
                gene_short_name = m.group(1)
                continue

            '''Translated protein product name attribute in CDS information'''
            m = re.match(PRODUCT_PATTERN, line)
            if m and is_in_cds_region:
                cds_region.product = m.group(1)

'''
Protein information parser class to extract protein information from the data handle returned from the
Entrez efetch function

---Attributes---
protein_info_set: Target protein information, a mapping between the target protein accession number and
                  its associated ProteinInfo object
'''
class ProteinXMLParser:
    '''Constructor'''
    def __init__(self):
        self._protein_info_set = dict()

    '''
    Function name: _parse_seq_definition
    Inputs       : Sequence definition data
    Outputs      : Protein name and source organism (can be empty string)
    Description  : Extracts the protein name and the source organism (an optional field) from the
                   sequence definition data
    '''
    @staticmethod
    def _parse_seq_definition(seq_def):
        m = re.match(r'^([^\[\]]+)(\s+\[(.+)\])?$', seq_def)
        if m:
            return m.group(1), m.group(3)
        else:
            return seq_def, ''

    '''
    Function name: parse
    Inputs       : Protein information data handle
    Outputs      : None
    Description  : Parses the data handle to extract target protein information. The data fetched is in
                   XML format
    '''
    def parse(self, handle):
        protein_info = ProteinInfo()

        xml_tree = ElementTree.parse(handle)
        root_node = xml_tree.getroot()

        for genbank_seq_node in root_node:
            '''Beginning of a new protein information'''
            protein_acc_num = genbank_seq_node.findtext('GBSeq_accession-version')
            if protein_acc_num is None:
                continue

            '''Save the last protein information'''
            if protein_info.acc_num is not None:
                self._protein_info_set[protein_info.acc_num] = protein_info

            protein_info = ProteinInfo()
            '''
            Protein accession number (or ID) and protein sequence (string) attributes in protein
            information
            '''
            protein_info.acc_num = protein_acc_num
            protein_info.seq_str = genbank_seq_node.findtext('GBSeq_sequence')

            '''Protein name and source organism attributes in protein information'''
            seq_def = genbank_seq_node.findtext('GBSeq_definition')
            protein_info.name, protein_info.organism = self._parse_seq_definition(seq_def)

            feature_table_node = genbank_seq_node.find('GBSeq_feature-table')
            if feature_table_node is None:
                continue

            for feature_node in feature_table_node:
                feature_key_value = feature_node.findtext('GBFeature_key')
                if feature_key_value == 'CDS':
                    '''Abbreviated gene name attribute in protein information'''
                    for feature_qual_node in feature_node.iter('GBQualifier'):
                        if feature_qual_node.findtext('GBQualifier_name') == 'gene':
                            protein_info.coding_gene_short_name = feature_qual_node.findtext('GBQualifier_value')
                            break
                elif feature_key_value == 'Protein':
                    '''Sequence completeness attribute in protein information'''
                    seq_range_value = feature_node.findtext('GBFeature_location')
                    if seq_range_value is not None:
                        m = re.match(PROTEIN_REGION_PATTERN, seq_range_value)
                        if m:
                            protein_info.is_5_partial = (m.group(1) == '<')
                            protein_info.is_3_partial = (m.group(3) == '>')

                    for feature_qual_node in feature_node.iter('GBQualifier'):
                        '''Protein name attribute in protein information'''
                        if feature_qual_node.findtext('GBQualifier_name') == 'product':
                            protein_info.name = feature_qual_node.findtext('GBQualifier_value')
                            break
                elif feature_key_value == 'source':
                    for feature_qual_node in feature_node.iter('GBQualifier'):
                        '''Source organism attribute in protein information'''
                        if feature_qual_node.findtext('GBQualifier_name') == 'organism':
                            protein_info.organism = feature_qual_node.findtext('GBQualifier_value')
                            break

        if protein_info.acc_num is not None:
            self._protein_info_set[protein_info.acc_num] = protein_info

    def get_protein_info_set(self):
        return self._protein_info_set

'''
Protein sequence parser class to extract protein sequences from the data handle returned from the
Entrez efetch function

---Attributes---
protein_seqs: Target protein sequences, a mapping between the target protein accession number and
              its associated sequence string
'''
class ProteinSeqParser:
    '''Constructor'''
    def __init__(self):
        self._protein_seqs = dict()

    def get_protein_seqs(self):
        return self._protein_seqs

    '''
    Function name: parse
    Inputs       : Protein sequence data handle
    Outputs      : None
    Description  : Parses the data handle to extract protein sequences
    '''
    def parse(self, handle):
        for seq_record in SeqIO.parse(handle, 'fasta'):
            protein_acc_num = extract_protein_acc_num(seq_record.description)
            self._protein_seqs[protein_acc_num] = seq_record.seq

'''
class ProteinSummaryParser:
    def __init__(self):
        self._protein_info = dict()

    def get_protein_info(self):
        return self._protein_info

    def parse(self, protein_docsum_handle):
        for protein_summary in Entrez.parse(protein_docsum_handle):
            protein_acc_num = protein_summary['AccessionVersion']
            self._protein_info[protein_acc_num] = ProteinInfo(protein_acc_num,
                                                              protein_summary['Title'],
                                                              protein_summary['Length'])
'''
