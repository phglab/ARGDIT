#!/usr/bin/env python3

from ArgditLib.ProcLog import ProcLog
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('seq_db_path', help = 'nucleotide/protein database FASTA file path')
parser.add_argument('replace_seq_file_path', help = 'FASTA file path for replacement sequences')
parser.add_argument('output_seq_db_path', help = 'output database file path')
args = parser.parse_args()

replace_seq_record_count = 0
with open(args.replace_seq_file_path, 'rU') as f:
    replace_seq_records = dict()
    for seq_record in SeqIO.parse(f, 'fasta'):
        replace_seq_records[seq_record.description] = seq_record
        replace_seq_record_count += 1

input_seq_record_count = 0
export_seq_record_count = 0
direct_export_count = 0
replace_count = 0

with open(args.output_seq_db_path, 'w') as fw:
    with open(args.seq_db_path, 'rU') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            input_seq_record_count += 1
            if seq_record.description in replace_seq_records:
                SeqIO.write(replace_seq_records[seq_record.description], fw, 'fasta')
                export_seq_record_count += 1
                replace_count += 1
            else:
                SeqIO.write(seq_record, fw, 'fasta')
                export_seq_record_count += 1
                direct_export_count += 1

summary = list()
summary_stmt = ProcLog.create_summary_stmt(input_seq_record_count, 'read from {}'.format(args.seq_db_path))
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(replace_seq_record_count, 'read from {}'.format(args.replace_seq_file_path),
                                           'replacement')
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(direct_export_count,
                                           'directly exported to {}'.format(args.output_seq_db_path), 'input')
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(replace_count, 'exported to {}'.format(args.output_seq_db_path),
                                           'replacement')
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(export_seq_record_count, 'exported in total')
summary.append(summary_stmt)
ProcLog.export_ext_summary(sys.stdout, summary)
