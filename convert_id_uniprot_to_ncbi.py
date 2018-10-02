#!/usr/bin/env python3

from ArgditLib.Config import Config
from ArgditLib.EntrezDBAccess import set_entrez_email, search_protein_seqs, search_protein_seqs_hist
from ArgditLib.ProcLog import ProcLog
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError, URLError
import argparse
import re
import sys
import urllib.parse
import urllib.request

UNIPROT_ACC_NUM_PATTERN = r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'
ID_CONVERT_URL = 'http://www.uniprot.org/uploadlists/'
MAX_ATTEMPT = 3

ProcLog.init_logs()

config = Config('config.ini')
if ProcLog.has_exec_error():
    ProcLog.export_exec_error(sys.stdout)
    sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument('seq_db_path', help = 'FASTA file path for protein database with Uniprot IDs')
parser.add_argument('output_seq_db_path',
                    help = 'output file path for protein database with converted NCBI protein accession no.')
args = parser.parse_args()

uniprot_ids = set()

with open(args.seq_db_path, 'rU') as f:
    for seq_record in SeqIO.parse(f, 'fasta'):
        m = re.search(UNIPROT_ACC_NUM_PATTERN, seq_record.description)
        if m:
            uniprot_ids.add(m.group(1))

params = {
    'from':'ACC',
    'to':'EMBL',
    'format':'tab',
    'query':' '.join(uniprot_ids)
}

print('Searching for Uniprot ID - NCBI accession number mapping...')

attempt = 0
while attempt < MAX_ATTEMPT:
    try:
        data = urllib.parse.urlencode(params).encode('utf-8')
        request = urllib.request.Request(ID_CONVERT_URL, data)
        request.add_header('User-Agent', 'Python %s' % config.entrez_email)
        response = urllib.request.urlopen(request)
        uniprot_to_ncbi_id_map_tbl = response.read().decode('utf-8')
        break
    except HTTPError as http_error:
        if 400 <= http_error.code <= 599:
            time.sleep(20)
            attempt += 1
            continue
        else:
            ProcLog.log_exec_error('UniProt: {}'.format(http_error.reason))
            break
    except URLError as url_error:
        ProcLog.log_exec_error('UniProt: {}'.format(str(url_error.reason)))
        break

if attempt == MAX_ATTEMPT:
    ProcLog.log_exec_error('Connecting UniProt 3 times and failed')

if ProcLog.has_exec_error():
    ProcLog.export_exec_error(sys.stdout)
    sys.exit()

uniprot_to_ncbi_ids_map = dict()
ncbi_protein_acc_nums = set()
uniprot_to_ncbi_id_map_pair_strs = uniprot_to_ncbi_id_map_tbl.splitlines()

for id_map_pair_str in uniprot_to_ncbi_id_map_pair_strs:
    id_map_pair = id_map_pair_str.split('\t')    
    if id_map_pair[0] in uniprot_to_ncbi_ids_map:
        uniprot_to_ncbi_ids_map[id_map_pair[0]].append(id_map_pair[1])
    else:
        uniprot_to_ncbi_ids_map[id_map_pair[0]] = [id_map_pair[1]]

    ncbi_protein_acc_nums.add(id_map_pair[1])

print('Loading data from NCBI...')

set_entrez_email(config.entrez_email)
if ProcLog.has_exec_error():
    ProcLog.export_exec_error(sys.stdout)
    sys.exit()

ncbi_protein_seqs = search_protein_seqs(ncbi_protein_acc_nums)

if ProcLog.has_exec_error():
    ProcLog.export_exec_error(sys.stdout)
    sys.exit()

input_seq_record_count = 0
export_seq_record_count = 0
no_uniprot_id_count = 0
unresolved_uniprot_id_count = 0

with open(args.output_seq_db_path, 'w') as fw:
    with open(args.seq_db_path, 'rU') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            input_seq_record_count += 1
            m = re.search(UNIPROT_ACC_NUM_PATTERN, seq_record.description)
            if m:
                uniprot_id = m.group(1)
            else:
                ProcLog.log_exec_msg('Uniprot ID not found in {}'.format(seq_record.description))
                continue

            if uniprot_id in uniprot_to_ncbi_ids_map:
                protein_seq_str = str(seq_record.seq)
                matched_ncbi_protein_acc_num = None
                for ncbi_protein_acc_num in uniprot_to_ncbi_ids_map[uniprot_id]:
                    if protein_seq_str == str(ncbi_protein_seqs[ncbi_protein_acc_num]):
                        matched_ncbi_protein_acc_num = ncbi_protein_acc_num
                        break

                if matched_ncbi_protein_acc_num is None:
                    matched_ncbi_protein_acc_num = uniprot_to_ncbi_ids_map[uniprot_id][0]

                ncbi_id_seq_header = seq_record.description.replace(uniprot_id, matched_ncbi_protein_acc_num)
                ncbi_id_seq_record = SeqRecord(seq_record.seq, id = ncbi_id_seq_header, name = '', description = '')
                SeqIO.write(ncbi_id_seq_record, fw, 'fasta')
                export_seq_record_count += 1
            else:
                ProcLog.log_exec_msg('NCBI accession no. not found for {}'.format(seq_record.description))
                unresolved_uniprot_id_count += 1

summary = list()
summary_stmt = ProcLog.create_summary_stmt(input_seq_record_count, 'read from {}'.format(args.seq_db_path))
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(no_uniprot_id_count, 'with no Uniprot ID')
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(unresolved_uniprot_id_count, 'with unresolved Uniprot ID')
summary.append(summary_stmt)
summary_stmt = ProcLog.create_summary_stmt(export_seq_record_count, 'exported to {}'.format(args.output_seq_db_path))
summary.append(summary_stmt)

ProcLog.export_exec_msg(sys.stdout, True)
ProcLog.export_ext_summary(sys.stdout, summary)
