#!/usr/bin/env python3

from ArgditLib.Utils import create_supp_file_path, check_seq_type
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os

def parse_seq_file(seq_db_path, diff_mode, is_check_rev_comp_seq):
    db_seq_records = dict()
    db_rev_seq_records = dict()
    num_of_seqs = 0

    with open(seq_db_path, 'rU') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            num_of_seqs += 1
            rev_seq_record_key = None

            if diff_mode == 'h':
                seq_record_key = seq_record.description
            elif diff_mode == 's':
                seq_record_key = str(seq_record.seq)
                if seq_record_key in db_rev_seq_records:
                    continue

                if is_check_rev_comp_seq and check_seq_type(seq_record_key) == 'NT':
                    rev_seq_record_key = str(seq_record.seq.reverse_complement())
            else:
                seq_str = str(seq_record.seq)
                seq_record_key = '{}###{}'.format(seq_record.description, seq_str)
                if seq_record_key in db_rev_seq_records:
                    continue

                if is_check_rev_comp_seq and check_seq_type(seq_str) == 'NT':
                    rev_seq_record_key = '{}###{}'.format(seq_record.description,
                                                          str(seq_record.seq.reverse_complement()))

            db_seq_records[seq_record_key] = seq_record
            if rev_seq_record_key is not None:
                db_rev_seq_records[rev_seq_record_key] = seq_record

    return db_seq_records, db_rev_seq_records, num_of_seqs

parser = argparse.ArgumentParser()
parser.add_argument('seq_db_path', help = 'nucleotide/protein database FASTA file path')
parser.add_argument('old_seq_db_path', help = 'FASTA file path for previous database release')
parser.add_argument('-m', '--mode', default = 'b', choices = ['h', 's', 'b'],
                    help = 'diff mode (h:header; s:sequence; b:both [default])')
parser.add_argument('-o', '--old', action = 'store_true', help = 'export also sequences found in old database only')
parser.add_argument('-f', '--forward', action = 'store_true',
                    help = 'do not check reverse complement nucleotide sequences')
args = parser.parse_args()

new_db_seq_records, new_db_rev_seq_records, num_of_new_db_seqs = parse_seq_file(args.seq_db_path, args.mode,
                                                                                not args.forward)
old_db_seq_records, old_db_rev_seq_records, num_of_old_db_seqs = parse_seq_file(args.old_seq_db_path, args.mode,
                                                                                not args.forward)

new_db_seq_keys = set(new_db_seq_records.keys())
old_db_seq_keys = set(old_db_seq_records.keys())
common_seq_keys = new_db_seq_keys & old_db_seq_keys

old_db_rev_seq_keys = set(old_db_rev_seq_records.keys())
common_seq_keys = ((new_db_seq_keys - common_seq_keys) & old_db_rev_seq_keys) | common_seq_keys

new_db_only_seq_records = list()
for new_db_only_seq_key in (new_db_seq_keys - common_seq_keys):
    new_db_only_seq_records.append(new_db_seq_records[new_db_only_seq_key])

new_db_only_seq_file_path = create_supp_file_path(args.seq_db_path, '_only.fa')
with open(new_db_only_seq_file_path, 'w') as f:
    SeqIO.write(new_db_only_seq_records, f, 'fasta')

common_seq_records = list()
for common_seq_key in common_seq_keys:
    common_seq_records.append(new_db_seq_records[common_seq_key])

common_seq_file_path = create_supp_file_path(args.seq_db_path, '_shared.fa')
with open(common_seq_file_path, 'w') as f:
    SeqIO.write(common_seq_records, f, 'fasta')

if args.old:
    common_seq_keys = new_db_seq_keys & old_db_seq_keys
    new_db_rev_seq_keys = set(new_db_rev_seq_records.keys())
    common_seq_keys = ((old_db_seq_keys - common_seq_keys) & new_db_rev_seq_keys) | common_seq_keys

    old_db_only_seq_records = list()
    for old_db_only_seq_key in (old_db_seq_keys - common_seq_keys):
        old_db_only_seq_records.append(old_db_seq_records[old_db_only_seq_key])

    old_db_only_seq_file_path = create_supp_file_path(args.old_seq_db_path, '_only.fa')
    with open(old_db_only_seq_file_path, 'w') as f:
        SeqIO.write(old_db_only_seq_records, f, 'fasta')

if args.mode == 'h':
    seq_db_summary_stmt = '  {} sequences read {}  {} of them are unique in terms of header'
elif args.mode == 's':
    seq_db_summary_stmt = '  {} sequences read {}  {} of them are unique in terms of sequence'
else:
    seq_db_summary_stmt = '  {} sequences read {}  {} of them are unique in terms of header and/or sequence'

print('----- Input summary -----')
print('For {}:'.format(args.seq_db_path))
print(seq_db_summary_stmt.format(num_of_new_db_seqs, os.linesep, len(new_db_seq_records)))
print('For {}:'.format(args.old_seq_db_path))
print(seq_db_summary_stmt.format(num_of_old_db_seqs, os.linesep, len(old_db_seq_records)))
print('----- Diff results -----')
print('{} sequences are common and exported to {}'.format(len(common_seq_records), common_seq_file_path))
print('{} sequences are present in {} only and exported to {}'.format(len(new_db_only_seq_records), args.seq_db_path,
                                                                      new_db_only_seq_file_path))

if args.old:
    print('{} sequences are present in {} only and exported to {}'.format(len(old_db_only_seq_records),
                                                                          args.old_seq_db_path,
                                                                          old_db_only_seq_file_path))
