import itertools
import logging
from collections import Counter
from SeqInfo import SeqInfo
from pamscore import PAM30_score


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

def get_S2_scores(msa_seq_strs, S2, score_func):
    msa_len = len(msa_seq_strs[0])
    S2_seqs = [msa_seq_strs[s_idx] for s_idx in sorted(S2)]
    S2Comp = set(range(len(msa_seq_strs))) - S2
    S2Comp_seqs = [msa_seq_strs[s_idx] for s_idx in sorted(S2Comp)]
    
    S2_scores = []
    curr = 0
    for col_idx in range(msa_len):
        curr += min([score_func(s1, s2, col_idx) for i, s1 in enumerate(S2_seqs) for s2 in S2_seqs[i:]])
        curr -= 2 * max([score_func(s, c, col_idx) for s in S2_seqs for c in S2Comp_seqs])
        curr = max(0, curr)
        S2_scores.append(curr)
    return S2_scores

def get_score_and_pos(msa_seq_strs, S2, score_func):
    S2_scores = get_S2_scores(msa_seq_strs, S2, score_func)
    max_score = max(S2_scores)
    max_score_end = S2_scores.index(max_score) + 1
    max_score_start = max_score_end - 1
    while max_score_start > 0 and S2_scores[max_score_start - 1] > 0:
        max_score_start -= 1
    return max_score, (max_score_start, max_score_end)


