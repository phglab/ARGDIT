'''Module for sequence translation functions used by ARGDIT tools'''

from Bio.Data import CodonTable
from Bio.Seq import Seq
from .CDSSeqSegment import CDSSeqSegment
from .CDSRegion import CDSRegion
from .ProcLog import ProcLog
from .ProteinInfo import ProteinInfo
from . import Utils

class Translate:
    '''
    Function name: generate_partial_codon_translate_map
    Inputs       : Translation table
    Outputs      : Translation map for the partial codons, if any
    Description  : For a particular genetic code (in terms of translation table here), search for any
                   possible partial condons that can be unambiguously translated into amino acid
    '''
    @staticmethod
    def generate_partial_codon_translate_map(translate_table):
        partial_codon_translate_map = dict()
        all_nucleotides = translate_table.nucleotide_alphabet

        for codon_nucleotide1 in all_nucleotides:
            for codon_nucleotide2 in all_nucleotides:
                partial_codon = '{}{}'.format(codon_nucleotide1, codon_nucleotide2)
                last_translated_aa = None
                is_same_translation = True

                for codon_nucleotide3 in all_nucleotides:
                    codon = '{}{}'.format(partial_codon, codon_nucleotide3)
                    codon_seq = Seq(codon)
                    translated_aa = str(codon_seq.translate(table = translate_table))

                    if last_translated_aa is None:
                        last_translated_aa = translated_aa
                    elif not last_translated_aa == translated_aa:
                        is_same_translation = False
                        break

                if is_same_translation and last_translated_aa is not None:
                    partial_codon_translate_map[partial_codon] = last_translated_aa

        return partial_codon_translate_map

    '''
    Function name: init
    Inputs       : Genetic code
    Outputs      : Nil
    Description  : Obtain the corresponding start and stop codons, and the translation map for the partial
                   condons according to the genetic code provided
    '''
    @classmethod
    def init(cls, genetic_code):
        try:
            translate_table = CodonTable.unambiguous_dna_by_id[genetic_code]
        except KeyError as e:
            ProcLog.log_exec_error('The genetic code {} is not valid'.format(genetic_code))
            return

        cls._genetic_code = genetic_code
        cls._start_codons = set(translate_table.start_codons)
        cls._stop_codons = set(translate_table.stop_codons)
        cls._partial_codon_translate_map = cls.generate_partial_codon_translate_map(translate_table)

    @classmethod
    def get_start_codons(cls):
        return cls._start_codons

    @classmethod
    def get_stop_codons(cls):
        return cls._stop_codons

    '''
    Function name: _translate_cds_seq_segment
    Inputs       : Candidate CDS sequence segment and its matched target CDS regions
    Outputs      : All possible translated protein sequences
    Description  : The candidate CDS sequence segment is translated according to its matched target CDS region.
                   Multiple protein sequences are sometimes returned since some partial codons maybe translated
    '''
    @classmethod
    def _translate_cds_seq_segment(cls, cds_seq_segment, cds_region):
        if cds_seq_segment.is_5_partial:
            translate_start = cds_region.codon_start - 1
        else:
            translate_start = 0

        translate_end = translate_start + (cds_seq_segment.length - translate_start) // 3 * 3

        cds_seq = Seq(cds_seq_segment.seq_str[translate_start:translate_end])

        protein_seq_str = str(cds_seq.translate(table = cls._genetic_code, to_stop = False))
        if (not cds_seq_segment.is_5_partial and cds_seq_segment.seq_str[translate_start:3] in cls._start_codons) or \
            cds_region.is_enforce_start_codon:
            protein_seq_str = '{}{}'.format('M', protein_seq_str[1:])

        if not cds_seq_segment.is_3_partial and protein_seq_str[-1] == '*':
            protein_seq_str = protein_seq_str[:-1]

        '''
        After removal of the * character at the end of the protein sequence, replace other * in the sequence with X
        This is to handle possble stop codon suppression, e.g. the intI2 gene in HM043577.1
        '''
        protein_seq_str = protein_seq_str.replace('*', 'X')

        if cds_seq_segment.is_3_partial and (len(protein_seq_str) * 3) == (translate_end - translate_start) and \
            (cds_seq_segment.length - translate_end) == 2:
            residue_partial_codon_str = cds_seq_segment.seq_str[-2:]
            if residue_partial_codon_str in cls._partial_codon_translate_map:
                ext_protein_seq_str = '{}{}'.format(protein_seq_str,
                                      cls._partial_codon_translate_map[residue_partial_codon_str])
                return {protein_seq_str, ext_protein_seq_str}

        return {protein_seq_str}

    '''
    Function name: search_correct_cds_translation
    Inputs       : Candidate CDS sequence segments, target CDS regions, protein information downloaded from NCBI
                   protein database, protein IDs in non-accession number format
    Outputs      : Target protein information, matched CDS sequence segment, matched CDS region, and the
                   non-accession number format protein ID to protein accession number mapping
    Description  : For every candidate sequence segments predicted, when it is matched with a target CDS region,
                   translation is performed with it accordingly. If the translated protein matches with any of
                   the target protein, the correct CDS sequence and region are then found and returned. If the
                   protein ID in the matched target CDS region is not in accession number format, it is mapped
                   back to the protein accession number
    '''
    @classmethod
    def search_correct_cds_translation(cls, candidate_cds_seq_segments, target_cds_regions, genbank_protein_info_set,
                                       non_acc_fmt_protein_ids):
        for candidate_cds_seq_segment in candidate_cds_seq_segments:
            for target_cds_region in target_cds_regions:
                if target_cds_region.protein_id not in genbank_protein_info_set and len(non_acc_fmt_protein_ids) == 0:
                    continue

                if (candidate_cds_seq_segment.length != target_cds_region.length) or \
                   (candidate_cds_seq_segment.is_3_partial is not target_cds_region.is_3_partial):
                    continue

                if (candidate_cds_seq_segment.is_5_partial is not target_cds_region.is_5_partial):
                    if not (candidate_cds_seq_segment.is_5_partial and target_cds_region.is_enforce_start_codon):
                        continue

                candidate_protein_seq_strs = cls._translate_cds_seq_segment(candidate_cds_seq_segment,
                                                                            target_cds_region)

                if target_cds_region.protein_id in genbank_protein_info_set:
                    target_protein_info = genbank_protein_info_set[target_cds_region.protein_id]
                    if target_protein_info.seq_str in candidate_protein_seq_strs:
                        return target_protein_info, candidate_cds_seq_segment, target_cds_region, None
                else:
                    matched_protein_id = Utils.match_non_acc_fmt_genbank_protein_seqs(candidate_protein_seq_strs,
                                                                                      genbank_protein_info_set,
                                                                                      non_acc_fmt_protein_ids)
                    if matched_protein_id is not None:
                        return genbank_protein_info_set[matched_protein_id], candidate_cds_seq_segment, \
                            target_cds_region, (target_cds_region.protein_id, matched_protein_id)

        return None, None, None, None

    '''
    Function name: _translate_gene_seq
    Inputs       : Gene sequence string
    Outputs      : Translated protein sequences
    Description  : Performs gene sequence translation for the entire input sequence (frame)
    '''
    @classmethod
    def _translate_gene_seq(cls, gene_seq):
        protein_seq_str = str(gene_seq.translate(table = cls._genetic_code, to_stop = True))
        translate_outputs = [protein_seq_str]

        if str(gene_seq[0:3]) in cls._start_codons:
            translate_outputs.append('M{}'.format(protein_seq_str[1:]))

        return translate_outputs

    '''
    Function name: _six_frame_translate
    Inputs       : Gene sequence string
    Outputs      : All translation products of six frames
    Description  : Performs six-frame translation, including partial codons when necessary. Translation
                   products of all the six frames are returned
    '''
    @classmethod
    def _six_frame_translate(cls, gene_seq_str):
        frame_outputs = list()
        extended_frame_outputs = list()
        full_seq_len = len(gene_seq_str)
        translate_lens = [full_seq_len // 3 * 3, (full_seq_len - 1) // 3 * 3, (full_seq_len - 2) // 3 * 3]

        forward_seq = Seq(gene_seq_str)
        rev_comp_seq = forward_seq.reverse_complement()

        for i in range(3):
            forward_translate_outputs = cls._translate_gene_seq(forward_seq[i:i + translate_lens[i]])
            rev_comp_translate_outputs = cls._translate_gene_seq(rev_comp_seq[i:i + translate_lens[i]])

            frame_outputs.append(forward_translate_outputs)
            frame_outputs.append(rev_comp_translate_outputs)

            if full_seq_len - (i + translate_lens[i]) == 2:
                residual_partial_codon_str = str(forward_seq[i + translate_lens[i]:full_seq_len])
                if (len(forward_translate_outputs[0]) * 3) == translate_lens[i] and \
                   residual_partial_codon_str in cls._partial_codon_translate_map:
                    extended_outputs = list()
                    for translate_output in forward_translate_outputs:
                        last_aa = cls._partial_codon_translate_map[residual_partial_codon_str]
                        extended_outputs.append('{}{}'.format(translate_output, last_aa))

                    extended_frame_outputs.append(extended_outputs)
                else:
                    extended_frame_outputs.append(None)

                residual_partial_codon_str = str(rev_comp_seq[i + translate_lens[i]:full_seq_len])
                if (len(rev_comp_translate_outputs[0]) * 3) == translate_lens[i] and \
                   residual_partial_codon_str in cls._partial_codon_translate_map:
                    extended_outputs = list()
                    for translate_output in rev_comp_translate_outputs:
                        last_aa = cls._partial_codon_translate_map[residual_partial_codon_str]
                        extended_outputs.append('{}{}'.format(translate_output, last_aa))

                    extended_frame_outputs.append(extended_outputs)
                else:
                    extended_frame_outputs.append(None)
            else:
                for j in range(2):
                    extended_frame_outputs.append(None)

        return frame_outputs, extended_frame_outputs

    '''
    Function name: translate
    Inputs       : Gene sequence string
    Outputs      : All frame translation prodicts from six-frame translation
    Description  : Captures and returns all possible translation products from each frame in six-frame
                   translation
    '''
    @classmethod
    def translate(cls, gene_seq_str):
        if len(gene_seq_str) < 3:
            return set()

        all_frame_outputs = set()
        frame_outputs, extended_frame_outputs = cls._six_frame_translate(gene_seq_str)

        for i in range(len(frame_outputs)):
            all_frame_outputs.update(frame_outputs[i])

            if extended_frame_outputs[i] is not None:
                all_frame_outputs.update(extended_frame_outputs[i])

        return all_frame_outputs
