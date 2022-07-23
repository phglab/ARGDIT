'''Multiple sequence alignment functions'''

from Bio import SeqIO
import os
import shlex
import subprocess
import tempfile

'''
Function name: multiSeqAlign
Inputs       : Protein sequences in FASTA format (as list of Bio.SeqRecord)
Outputs      : Aligned sequences in FASTA format (as list of Bio.SeqRecord)
Description  : Performs multiple sequence alignment using MUSCLE
'''
def muscleAlign(seq_records, align_file_path = None):
    with tempfile.NamedTemporaryFile(mode = 'w', delete = False) as f:
        SeqIO.write(seq_records, f, 'fasta')
        tmp_seq_file_path = f.name

    if align_file_path is None:
        with tempfile.NamedTemporaryFile(mode = 'w', delete = False) as f:
            tmp_align_file_path = f.name
    else:
        tmp_align_file_path = align_file_path

    muscle_cmd_line = f'muscle -align {tmp_seq_file_path} -output {tmp_align_file_path}'

    subprocess.run(shlex.split(muscle_cmd_line),
                   stdin = subprocess.DEVNULL,
                   stdout = subprocess.DEVNULL,
                   stderr = subprocess.DEVNULL,
                   universal_newlines=True)

    os.remove(tmp_seq_file_path)

    return tmp_align_file_path
