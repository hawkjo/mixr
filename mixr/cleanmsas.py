import os
import itertools
import logging
from collections import Counter
from Bio import SeqIO
from SeqInfo import SeqInfo
from pamscore import build_PAM30_score_func


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

def python_intervals_overlap(interval1, interval2):
    query_start, query_end = interval1
    subject_start, subject_end = interval2
    return (query_start <= subject_start < query_end
            or query_start < subject_end <= query_end
            or subject_start <= query_start < subject_end
            or subject_start < query_end <= subject_end)

def keep_only_good_exons(rec, good_exons):
    start, end = good_exons[0]
    new_rec = rec[start:end]
    for start, end in good_exons[1:]:
        new_rec += rec[start:end]
    return new_rec

def update_exons(exons, good_exons):
    last_end = 0
    to_subtract = 0
    new_exons = []
    for start, end in good_exons:
        to_subtract += start - last_end
        new_exons.append((start - to_subtract, end - to_subtract))
        last_end = end
    return new_exons

def next_significant_score(msa_recs,
                           seq_info,
                           score_func,
                           shortest_allowed,
                           largest_to_smallest_subset=True):

    species_names = [str(rec.id) for rec in msa_recs]
    msa_seq_strs = [str(rec.seq) for rec in msa_recs]
    msa_len = len(msa_seq_strs[0])
    neg_ctrl_len = 100 * msa_len
    loc_S2s = get_second_consensus_subsets(msa_seq_strs)

    for S2 in sorted(loc_S2s, key=len, reverse=largest_to_smallest_subset):  
        score, pos_tup = get_score_and_pos(msa_seq_strs, S2, score_func)
        if pos_tup[1] - pos_tup[0] < shortest_allowed:
            continue

        neg_ctrl_seq_strs = seq_info.get_negative_control_seqs(species_names, neg_ctrl_len)
        neg_ctrl_score, neg_ctrl_pos_tup = get_score_and_pos(neg_ctrl_seq_strs, S2, score_func)

        if score > neg_ctrl_score:
            return score, neg_ctrl_score, pos_tup, S2
    return None

def filter_exons_then_species(fname,
                              seq_info,
                              score_func,
                              shortest_allowed=3):

    # Filter sequences, attempting to keep as many species as possible. Specifically:
    #
    # while still badness:
    #     store pre-exon-removal state
    #     while MSA remaining and still badness:
    #         remove largest significant-score subset exons
    #     if whole MSA removed:
    #         restore pre-exon-removal state
    #         remove smallest significant-score subset species
    # return results
    
    msa_recs = seq_info.msa_recs_given_fname[fname]
    cds_msa_recs = seq_info.cds_msa_recs_given_fname[fname]
    exons = seq_info.exons_given_fname[fname]
    removal_actions = []
    
    still_badness = True
    while still_badness:
        pre_exon_removal_msa_recs = msa_recs
        pre_exon_removal_cds_msa_recs = cds_msa_recs
        pre_exon_removal_removal_actions = removal_actions[:]
        pre_exon_removal_exons = exons[:]
        whole_MSA_removed = False
        
        while True:
            res = next_significant_score(msa_recs,
                                         seq_info,
                                         score_func,
                                         shortest_allowed,
                                         largest_to_smallest_subset=True)
            if res is None:
                still_badness = False
                break
                
            score, neg_ctrl_score, pos_tup, S2 = res
            species_names = [str(rec.id) for rec in msa_recs]
            bad_species_names = [species_names[i] for i in S2]
            
            good_exons = [exon for exon in exons if not python_intervals_overlap(pos_tup, exon)]
            if len(good_exons) == 0:
                whole_MSA_removed = True
                break
                
            bad_exons = [exon for exon in exons if exon not in good_exons]
            exons = update_exons(exons, good_exons)
            
            msa_recs = [keep_only_good_exons(rec, good_exons) for rec in msa_recs]
            cds_good_exons = [(3 * start, 3 * end) for start, end in good_exons]
            cds_msa_recs = [keep_only_good_exons(rec, cds_good_exons)
                            for rec in cds_msa_recs]

            
            removal_actions.append(
                ['exon removal',
                 score,
                 neg_ctrl_score,
                 pos_tup,
                 bad_species_names,
                 bad_exons,
                 good_exons,
                 exons]
            )

        if whole_MSA_removed:
            msa_recs = pre_exon_removal_msa_recs
            cds_msa_recs = pre_exon_removal_cds_msa_recs
            removal_actions = pre_exon_removal_removal_actions[:]
            exons = pre_exon_removal_exons[:]
            
            res = next_significant_score(msa_recs,
                                         seq_info,
                                         score_func,
                                         shortest_allowed,
                                         largest_to_smallest_subset=False)
            if res is None:
                break

            score, neg_ctrl_score, pos_tup, S2 = res
            species_names = [str(rec.id) for rec in msa_recs]
            bad_species_names = [species_names[i] for i in S2]
            msa_recs = [rec for rec in msa_recs if str(rec.id) not in bad_species_names]
            cds_msa_recs = [rec for rec in cds_msa_recs if str(rec.id) not in bad_species_names]
            removal_actions.append(
                ['species removal',
                 score,
                 neg_ctrl_score,
                 pos_tup,
                 bad_species_names]
            )
                
    return msa_recs, cds_msa_recs, removal_actions

def remove_gaps(fname, msa_recs, cds_msa_recs, prefilter_exons, removal_actions):
    # Find post-filter exons
    exons = prefilter_exons
    for action in removal_actions:
        if action[0].startswith('exon'):
            exons = action[-1]

    # Find all-gap columns
    msa_seq_strs = [str(rec.seq) for rec in msa_recs]
    cds_msa_seq_strs = [str(rec.seq) for rec in cds_msa_recs]
    msa_len = len(msa_seq_strs[0])
    all_gap_col_idxs = []
    for col_idx in range(msa_len):
        if all([seq[col_idx] == '-' for seq in msa_seq_strs]):
            for cds_col_idx in range(3 * col_idx, 3 * col_idx + 3):
                assert all([seq[cds_col_idx] == '-' for seq in cds_msa_seq_strs]), fname
            all_gap_col_idxs.append(col_idx)
            
    # Remove gaps
    for col_idx in sorted(all_gap_col_idxs, reverse=True):  # Remove gaps right to left so idxs work
        msa_recs = [rec[:col_idx] + rec[col_idx+1:] for rec in msa_recs]
        new_len = len(msa_recs[0])
        cds_col_idx = 3 * col_idx
        cds_msa_recs = [rec[:cds_col_idx] + rec[cds_col_idx+3:] for rec in cds_msa_recs]
        new_exons = []
        for start, end in exons:
            if col_idx < start:
                start -= 1
            if col_idx < end:
                end -= 1
            new_exons.append((start, end))
        exons = new_exons
        
    # Check that everything looks good if things changed
    if all_gap_col_idxs:
        for rec, cds_rec in zip(msa_recs, cds_msa_recs):
            assert 3 * len(rec) == len(cds_rec), fname
        
        assert len(set(map(len, msa_recs))) == 1, fname
        new_msa_len = len(msa_recs[0])
        assert new_msa_len == msa_len - len(all_gap_col_idxs), (fname, msa_len, new_msa_len, all_gap_col_idxs, msa_recs)
        assert new_msa_len == exons[-1][-1], (fname, msa_len, new_msa_len, all_gap_col_idxs, exons)
        for (start1, end1), (start2, end2) in zip(exons[:-1], exons[1:]):
            assert end1 == start2, fname
        
        for col_idx in range(new_msa_len):
            msa_seq_strs = [str(rec.seq) for rec in msa_recs]
            assert not all([seq[col_idx] == '-' for seq in msa_seq_strs]), fname
            
    return msa_recs, cds_msa_recs, exons, all_gap_col_idxs

def exon_pos_from_exons(exons):
    return [start for start, end in exons] + [exons[-1][-1]]

def clean_msas(arguments):
    seq_info = SeqInfo(arguments)
    PAM30_score = build_PAM30_score_func()

    for out_dir in [arguments.out_prot_dir, arguments.out_cds_dir]:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    out_exon_pos_given_fname = {}
    for fname in seq_info.fnames:
        log.info('Cleaning {}'.format(fname))
        msa_recs, cds_msa_recs, removal_actions = filter_exons_then_species(
            fname,
            seq_info,
            PAM30_score
        )
        msa_recs, cds_msa_recs, exons, all_gap_col_idxs = remove_gaps(
            fname,
            msa_recs,
            cds_msa_recs,
            seq_info.exons_given_fname[fname],
            removal_actions
        )
        if not removal_actions:
            log.info('\tAlready clean'.format(fname))
        else:
            for actions in removal_actions:
                if actions[0].startswith('exon'):
                    log.info('\t{}: {}'.format(actions[0], str(actions[5])))
                else:
                    assert actions[0].startswith('species'), actions
                    log.info('\t{}: {}'.format(actions[0], str(actions[-1])))
        
        SeqIO.write(msa_recs, open(os.path.join(arguments.out_prot_dir, fname), 'w'), 'fasta')
        SeqIO.write(cds_msa_recs, open(os.path.join(arguments.out_cds_dir, fname), 'w'), 'fasta')
        out_exon_pos_given_fname[fname] = exon_pos_from_exons(exons)

    with open(arguments.out_exon_pos_file, 'w') as out:
        for fname, exon_pos in sorted(out_exon_pos_given_fname.items()):
            out.write('{}\t{}\n'.format(fname, '\t'.join(map(str, exon_pos))))

