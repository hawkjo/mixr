import itertools
import logging
from collections import Counter
from SeqInfo import SeqInfo


log = logging.getLogger(__name__)

def get_second_consensus_subsets(msa_seq_strs):
    msa_len = len(msa_seq_strs[0])
    second_consensus_subsets = set()
    for col_idx in range(msa_len):
        col_cntr = Counter([seq[col_idx] for seq in msa_seq_strs])
        col_subsets = set()
        for base, count in col_cntr.items():
            if count <= float(len(msa_seq_strs)) / 2.0:
                col_subsets.add(
                    frozenset([j for j, seq in enumerate(msa_seq_strs) if seq[col_idx] == base])
                )
        if not col_subsets:
            continue
        second_consensus_subsets.update(col_subsets)
        for n in range(2, len(col_subsets) + 1):
            for subs in itertools.combinations(col_subsets, n):
                new_sub = frozenset(el for sub in subs for el in sub)
                if len(new_sub) <= float(len(msa_seq_strs)) / 2.0:
                    second_consensus_subsets.add(new_sub)
    return second_consensus_subsets
