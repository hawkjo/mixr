import os
import logging
from Bio import SeqIO

log = logging.getLogger(__name__)

class SeqInfo(object):
    """
    A container class for MSA fnames, records, and exon info.
    """
    def __init__(self, arguments):
        self.arguments = arguments
        self.process_exon_pos_file()
        self.load_msas()
        self.remove_singleton_msas()

    def process_exon_pos_file(self):
        exon_pos_given_fname = {}
        for line in open(self.arguments.exon_pos_fpath):
            words = line.strip().split()
            exon_pos_given_fname[words[0]] = map(int, words[1:])

        self.exons_given_fname = {fname: zip(exon_pos[:-1], exon_pos[1:])
                                  for fname, exon_pos in exon_pos_given_fname.items()}

        self.fnames = self.exons_given_fname.keys()
        self.fnames.sort()

    def load_msas(self):
        self.msa_recs_given_fname = {}
        for fname in self.fnames:
            fpath = os.path.join(arguments.in_dir, fname)
            self.msa_recs_given_fname[fname] = list(SeqIO.parse(open(fpath), 'fasta'))

    def remove_singleton_msas(self):
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

