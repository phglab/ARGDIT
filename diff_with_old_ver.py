#!/usr/bin/env python3

from ArgditLib.Utils import create_supp_file_path
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('seq_db_path', help = 'nucleotide/protein database FASTA file path')
parser.add_argument('old_seq_db_path', help = 'FASTA file path for previous database release')
args = parser.parse_args()

new_db_seq_records = dict()
old_db_seq_records = dict()

with open(args.seq_db_path, 'rU') as f:
    for seq_record in SeqIO.parse(f, 'fasta'):
        new_db_seq_records[seq_record.description] = seq_record

with open(args.old_seq_db_path, 'rU') as f:
    for seq_record in SeqIO.parse(f, 'fasta'):
        old_db_seq_records[seq_record.description] = seq_record

new_db_seq_hdrs = set(new_db_seq_records.keys())
old_db_seq_hdrs = set(old_db_seq_records.keys())
common_seq_hdrs = new_db_seq_hdrs & old_db_seq_hdrs

new_db_only_seq_records = list()
for new_db_only_seq_hdr in (new_db_seq_hdrs - common_seq_hdrs):
    new_db_only_seq_records.append(new_db_seq_records[new_db_only_seq_hdr])

common_seq_records = list()
for common_seq_hdr in common_seq_hdrs:
    new_db_seq_record= new_db_seq_records[common_seq_hdr]
    old_db_seq_record = old_db_seq_records[common_seq_hdr]
    if str(new_db_seq_record.seq) == str(old_db_seq_record.seq):
        common_seq_records.append(new_db_seq_record)
    else:
        new_db_only_seq_records.append(new_db_seq_record)

new_db_only_seq_file_path = create_supp_file_path(args.seq_db_path, '_only.fa')
with open(new_db_only_seq_file_path, 'w') as f:
    SeqIO.write(new_db_only_seq_records, f, 'fasta')

common_seq_file_path = create_supp_file_path(args.seq_db_path, '_shared.fa')
with open(common_seq_file_path, 'w') as f:
    SeqIO.write(common_seq_records, f, 'fasta')

print('{} sequences read from {}'.format(len(new_db_seq_records), args.seq_db_path))
print('{} sequences read from {}'.format(len(old_db_seq_records), args.old_seq_db_path))
print('{} sequences are present in {} only and exported to {}'.format(len(new_db_only_seq_records), args.seq_db_path,
                                                                      new_db_only_seq_file_path))
print('{} sequences are shared with {} and exported to {}'.format(len(common_seq_records), args.old_seq_db_path,
                                                                  common_seq_file_path))
