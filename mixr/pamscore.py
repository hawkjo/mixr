import os
import logging
from collections import Counter


log = logging.getLogger(__name__)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

pam_fpath = os.path.join(os.path.abspath(os.path.join(THIS_DIR, '..')), 'resources', 'pam30.txt')

log.info('Loading PAM30 scoring matrix')
with open(pam_fpath) as f:
    line = next(f).strip()
    while line.startswith('#'):
        line = next(f).strip()
    aas = line.split()
    PAM30 = {aa: {} for aa in aas}
    for line in f:
        words = line.strip().split()
        from_aa = words[0]
        for to_aa, score_str in zip(aas, words[1:]):
            PAM30[from_aa][to_aa] = int(score_str)
assert all(PAM30[aa1][aa2] == PAM30[aa2][aa1] for aa1 in aas for aa2 in aas)

all_PAM_scores = []
for i, from_aa in enumerate(aas[:-5]):
    for to_aa in aas[i:-5]:
        all_PAM_scores.append(PAM30[from_aa][to_aa])
        
PAM30_gap_open_score = max(Counter(all_PAM_scores))
PAM30_gap_extension_score = -1
            
def PAM30_score(seq1, seq2, col_idx):
    if seq1[col_idx] == '-' == seq2[col_idx]:
        return 0
    for seq in [seq1, seq2]:
        if seq[col_idx] is None:
            return PAM30_gap_open_score
        if seq[col_idx] == '-':
            if col_idx > 0 and seq[col_idx - 1] == '-':
                return PAM30_gap_extension_score
            else:
                return PAM30_gap_open_score
    return PAM30[seq1[col_idx]][seq2[col_idx]]
