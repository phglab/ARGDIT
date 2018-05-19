'''Multiple sequence alignment functions'''

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
import subprocess
import sys

'''
Function name: multiSeqAlign
Inputs       : Protein sequences in FASTA format (as list of Bio.SeqRecord)
Outputs      : Aligned sequences in FASTA format (as list of Bio.SeqRecord)
Description  : Performs multiple sequence alignment using MUSCLE
'''
def muscleAlign(seq_records):
    muscle_cmd_line = MuscleCommandline()

    child_process = subprocess.Popen(str(muscle_cmd_line),
                                     stdin = subprocess.PIPE,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.DEVNULL,
                                     universal_newlines = True)
    SeqIO.write(seq_records, child_process.stdin, 'fasta')
    child_process.stdin.close()
    
    return AlignIO.read(child_process.stdout, 'fasta')
