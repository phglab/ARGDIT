from .MultiSeqAlign import muscleAlign
from Bio import AlignIO
from Bio import SeqIO
from multiprocessing import Pool
import os
import re
import subprocess
import tempfile

'''Ontology annotator class to perform ontology annotation for ARG protein sequences'''
class OntologyAnnotator:
    '''
    Function name: create_hmm
    Inputs       : 2-tuple containing an ontology label and its underlying protein sequence records
    Outputs      : File path of the created HMM file
    Description  : First performs multiple sequence alignment of the input protein sequences, and then
                   builds the HMM from the alignment results. The HMM is generated as a file in the
                   temporary directory
    '''
    @staticmethod
    def create_hmm(otl_label_protein_seq_record_tuple):
        aligned_protein_seq_records = muscleAlign(otl_label_protein_seq_record_tuple[1])
        protein_seq_align_file_path = os.path.join(tempfile.gettempdir(), otl_label_protein_seq_record_tuple[0])
        
        with open(protein_seq_align_file_path, 'w') as f:
            AlignIO.write(aligned_protein_seq_records, f, 'fasta')

        with tempfile.NamedTemporaryFile(delete = False) as ftemp:
            hmm_file_path = ftemp.name

        subprocess_args = ['hmmbuild', hmm_file_path, protein_seq_align_file_path]
        child_process = subprocess.Popen(subprocess_args,
                                         stdout = subprocess.DEVNULL,
                                         stderr = subprocess.DEVNULL,
                                         universal_newlines = True)

        while (child_process.poll() is None):
            continue

        if os.path.exists(protein_seq_align_file_path):
            os.remove(protein_seq_align_file_path)

        return hmm_file_path

    '''
    Function name: create_compressed_hmm_db
    Inputs       : HMM file paths
    Outputs      : File path of the compressed HMM database
    Description  : Compresses all HMMs into a single HMM database
    '''
    @staticmethod
    def create_compressed_hmm_db(hmm_file_paths):
        with tempfile.NamedTemporaryFile(mode = 'wb', delete = False) as ftemp:
            hmm_db_file_path = ftemp.name

            for hmm_file_path in hmm_file_paths:
                with open(hmm_file_path, 'rbU') as f:
                    ftemp.write(f.read())

        subprocess_args = ['hmmpress', '-f', hmm_db_file_path]
        child_process = subprocess.Popen(subprocess_args,
                                         stdout = subprocess.DEVNULL,
                                         stderr = subprocess.DEVNULL,
                                         universal_newlines = True)

        while (child_process.poll() is None):
            continue

        if os.path.exists(hmm_db_file_path):
            os.remove(hmm_db_file_path)

            for hmm_file_path in hmm_file_paths:
                os.remove(hmm_file_path)

        return hmm_db_file_path

    '''Constructor to create compressed HMM database from the ARG protein sequences'''
    def __init__(self, otl_annot_protein_seq_record_grps):
        cpu_count = len(os.sched_getaffinity(0))

        with Pool(cpu_count) as pool:
            hmm_file_paths = list(pool.imap_unordered(self.create_hmm, otl_annot_protein_seq_record_grps.items()))

        pool.join()
        
        self.hmm_db_file_path = self.create_compressed_hmm_db(hmm_file_paths)

    '''
    Function name: parse_hmm_scan_output
    Inputs       : hmmscan result file path
    Outputs      : A sequence header to ontology class label map
    Description  : Parses the hmmscan results of the target ARG protein sequences to return their
                   ontology class labels
    '''
    @staticmethod
    def parse_hmm_scan_output(hmm_scan_output_tbl_path):
        otl_class_label_map = dict()

        with open(hmm_scan_output_tbl_path, 'rU') as f:
            last_query_name = None

            for line in f.readlines():
                if line.startswith('#'):
                    continue

                m = re.match(r'^(\S+)\s+\S+\s+(\S+).+$', line)
                if m:
                    if last_query_name is None or last_query_name != m.group(2):
                        otl_class_label_map[m.group(2)] = m.group(1)
                        last_query_name = m.group(2)

        os.remove(hmm_scan_output_tbl_path)

        return otl_class_label_map

    '''
    Function name: annotate_ontology_label
    Inputs       : Protein sequence records
    Outputs      : A sequence header to ontology class label map
    Description  : Peforms hmmscan function to predicted the ontology class labels for the input protein
                   sequences, and then parses the scan result file to return the predicted labels
    '''
    def annotate_ontology_label(self, protein_seq_records):
        with tempfile.NamedTemporaryFile(delete = False) as ftemp:
            hmm_scan_output_tbl_path = ftemp.name

        cpu_count = len(os.sched_getaffinity(0))
        subprocess_args = ['hmmscan', '--tblout', hmm_scan_output_tbl_path, '--cpu', str(cpu_count),
                           self.hmm_db_file_path, '-']
        child_process = subprocess.Popen(subprocess_args,
                                         stdin = subprocess.PIPE,
                                         stdout = subprocess.DEVNULL,
                                         stderr = subprocess.DEVNULL,
                                         universal_newlines = True)

        SeqIO.write(protein_seq_records, child_process.stdin, 'fasta')
        child_process.stdin.close()

        while (child_process.poll() is None):
            continue

        return self.parse_hmm_scan_output(hmm_scan_output_tbl_path)

    '''Destructor to release resources'''
    def __del__(self):
        hmm_db_file_basename = os.path.basename(self.hmm_db_file_path)
        temp_dir_path = tempfile.gettempdir()
        
        for temp_file_name in os.listdir(temp_dir_path):
            temp_file_name_parts = temp_file_name.rpartition('.')
            if temp_file_name_parts[0] == hmm_db_file_basename:
                os.remove(os.path.join(temp_dir_path, temp_file_name))
