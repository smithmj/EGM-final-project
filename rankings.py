import heapq
import itertools
import numpy as np
import time_machine

def rank_seq(filename, k, initial_state, prob_model, eq = False, long_period = False, immigration_prob = 0):
    '''
    Finds the sequence of antibiotics of length k optimizing the probability of 
    returning to the susceptible genotype 
    '''
    ab_seqs = list(itertools.permutations(time_machine.ANTIBIOTICS, k))
    h = []
    heapq.heapify(h)
    susc_index = time_machine.GENO.index('0000')
    for seq in ab_seqs:
        seq_matrices = []
        for ab in seq:
            m = prob_model(ab, immingration_prob = immigration_prob)
            if long_period:
                m = time_machine.eq_cycle_prob(m)
            seq_matrices.append(m)
        M = time_machine.cycle_prob(seq_matrices)
        if eq:
            M = time_machine.eq_cycle_prob(M)
        prob = np.dot(initial_state, M)[susc_index]
        # use the negative probability in the heap because the heap puts the 
        # smallest items on top. This means that we won't have to sort the 
        # heap later
        heapq.heappush(h, (-prob, seq))

    f = open(filename, 'w')
    f.write('Rank\t Probability\t Sequence\n')
    rank = 1
    while h != []:
        prob, seq = heapq.heappop(h)
        f.write('{}\t {}\t {}\n'.format(rank, -prob, seq))
        rank += 1
    f.close()
    return

def seq_probs(filename, k, initial_state, prob_model, eq = False, long_period = False, immigration_prob = 0):
    ab_seqs = list(itertools.permutations(time_machine.ANTIBIOTICS, k))
    susc_index = time_machine.GENO.index('0000')
    f = open(filename, 'w')
    f.write('Sequence\tProbability\n')
    for seq in ab_seqs:
        seq_matrices = []
        for ab in seq:
            m = prob_model(ab, immigration_prob = immigration_prob)
            if long_period:
                m = time_machine.eq_cycle_prob(m)
            seq_matrices.append(m)
        M = time_machine.cycle_prob(seq_matrices)
        if eq:
            M = time_machine.eq_cycle_prob(M)
        prob = np.dot(initial_state, M)[susc_index]
        f.write('{}\t{}\n'.format(seq, prob))
    f.close()
    return




if __name__ == "__main__":
    initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    
    seq_probs('probs_5perc_imm_2.txt', 2, initial_state, time_machine.cpm, immigration_prob = 0.05)
    seq_probs('probs_5perc_imm_eq_2.txt', 2, initial_state, time_machine.cpm, eq = True, immigration_prob = 0.05)
    seq_probs('probs_5perc_imm_long_period_2.txt', 2, initial_state, time_machine.cpm, long_period = True, immigration_prob = 0.05)

