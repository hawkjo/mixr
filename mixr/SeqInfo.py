import os
import random
import math
import logging
from Bio import SeqIO
from copy import deepcopy

log = logging.getLogger(__name__)

class SeqInfo(object):
    """
    A class for MSA fnames, records, exon info, and neg control seqs.
    """
    def __init__(self, arguments):
        self.arguments = arguments
        self.process_exon_pos_file()
        self.load_msas()
        self.remove_singleton_msas()
        self.infer_species()
        self.build_negative_control_seqs()

    def process_exon_pos_file(self):
        log.info('Reading exon pos file')
        exon_pos_given_fname = {}
        for line in open(self.arguments.exon_pos_fpath):
            words = line.strip().split()
            exon_pos_given_fname[words[0]] = map(int, words[1:])

        self.exons_given_fname = {fname: zip(exon_pos[:-1], exon_pos[1:])
                                  for fname, exon_pos in exon_pos_given_fname.items()}

        self.fnames = self.exons_given_fname.keys()
        self.fnames.sort()

    def load_msas(self):
        log.info('Loading MSAs')
        self.msa_recs_given_fname = {}
        for fname in self.fnames:
            fpath = os.path.join(arguments.in_dir, fname)
            self.msa_recs_given_fname[fname] = list(SeqIO.parse(open(fpath), 'fasta'))

    def remove_singleton_msas(self):
        log.info('Removing single-sequence MSAs')
        bad_fnames = set([fname for fname in self.fnames
                          if len(self.msa_recs_given_fname[fname]) < 2])
        self.fnames = [fname for fname in self.fnames if fname not in bad_fnames]
        for fname in bad_fnames:
            del self.msa_recs_given_fname[fname]
            del self.exons_given_fname[fname]


            bad_fnames = set([fname for fname in fnames if len(msa_recs_given_fname[fname]) < 2])
            fnames = [fname for fname in fnames if fname not in bad_fnames]
            for fname in bad_fnames:
                    del msa_recs_given_fname[fname]
                        del exon_pos_given_fname[fname]

    def infer_species(self):
        self.species = set()
        for msa_recs in self.msa_recs_given_fname.values():
            self.species.update([str(rec.id) for rec in msa_recs])
        self.sorted_species = list(self.species)
        self.sorted_species.sort()

    def build_negative_control_seqs(self):
        log.info('Building negative control sequences')

        all_species_cols = []
        for msa_recs in self.msa_recs_given_fname.values():
            if len(msa_recs) < len(self.species):
                continue
            msa_recs.sort(key=lambda rec: str(rec.id))
            msa_seq_strs = [str(rec.seq) for rec in msa_recs]
            msa_len = len(msa_seq_strs[0])
            for col_idx in range(msa_len):
                all_species_cols.append([seq[col_idx] for seq in msa_seq_strs])
        log.info('Number of all-species columns: {:,d}'.format(len(all_species_cols)))

        max_msa_len = max(len(msa_recs[0]) for msa_recs in self.msa_recs_given_fname.values())
        num_shuffles = int(math.ceil(float(100 * max_msa_len)/len(all_species_cols)))

        self.neg_control_seqs = ['' for _ in range(len(self.species))]
        for _ in range(num_shuffles):
            all_species_cols_copy = deepcopy(all_species_cols)
            random.shuffle(all_species_cols_copy)
            for row in range(len(self.neg_control_seqs)):
                self.neg_control_seqs[row] += ''.join([col[row] for col in all_species_cols_copy])
        assert len(set(map(len, self.neg_control_seqs))) == 1
        assert len(self.neg_control_seqs) == len(self.species)
        log.info('Negative control max sequence length: {:,d}   ({:.1f} x longest seq)'.format(
            len(self.neg_control_seqs[0]),
            len(self.neg_control_seqs[0])/float(max_msa_len)
        ))

    def get_negative_control_seqs(self, species_names, seqlen):
        return [seq[:seqlen] for species_name, seq in zip(self.sorted_species, neg_control_seqs)
                if species_name in species_names]
