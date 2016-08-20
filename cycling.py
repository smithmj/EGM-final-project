import itertools
import random
import numpy as np
import pandas as pd
import copy
import csv
import os
import re
import matplotlib.pyplot as plt  

ANTIBIOTICS = ['AMP','AM','CEC','CTX','ZOX','CXM','CRO','AMC', \
               'CAZ','CTT','SAM','CPR','CPD','TZP','FEP']

GENO = ['0000','1000','0100','0010',
        '0001','1100','1010','1001',
        '0110','0101','0011','1110',
        '1101','1011','0111','1111']

ADJ_MAT_DIR =  os.getcwd() + '/adj_mat'

GROWTH_RATE_FILE = os.getcwd() + '/growth_rate.csv'

def import_arr(filename):
    '''
    Imports a csv adjacency matrix file into a numpy array 
    '''
    with open(filename) as f:
        data = csv.reader(f)
        adj_mat = np.empty((0, 16))
        for row in data:
            row_arr = np.array([row]).astype(float)
            adj_mat = np.concatenate((adj_mat, row_arr))
    return adj_mat

def epm(ab, immigration_prob = 0):
    '''
    Creates the transition matrix for a given adjacency matrix 
    based on the equal probability model

    Default parameter immigration_prob will be added to each of the
    transition probabilities and each of the rows will be
    re-normalized
    '''
    adj_mat_file = '{}/{}_adj_mat.csv'.format(ADJ_MAT_DIR, ab)
    subt_mat = import_arr(adj_mat_file)
    for i in range(subt_mat.shape[0]):
        s = subt_mat[i].sum()
        if s != 0:
            subt_mat[i] = subt_mat[i] / s
    # add immigration probabilities
    if immigration_prob != 0:
        for i in range(subt_mat.shape[0]):
            num_non_zero = np.count_nonzero(subt_mat[i,:])
            num_zeros = subt_mat.shape[0] - num_non_zero
            for j in range(subt_mat.shape[1]):
                if subt_mat[i, j] == 0:
                    subt_mat[i, j] += immigration_prob / num_zeros
                else:
                    subt_mat[i, j] -= immigration_prob / num_non_zero
    return subt_mat

def cpm(ab, immigration_prob = 0):
    '''
    Creates the transition matrix for a given adjacency matrix 
    based on the correlated probability model
    '''
    adj_mat_file = '{}/{}_adj_mat.csv'.format(ADJ_MAT_DIR, ab)
    adj_mat = import_arr(adj_mat_file)
    ab = re.findall('[A-Z]{2,3}', adj_mat_file)[-1]
    growth_rate = import_arr(GROWTH_RATE_FILE)[ANTIBIOTICS.index(ab)]
    subt_mat = np.zeros(adj_mat.shape)
    for u in range(subt_mat.shape[0]):
        for v in range(subt_mat.shape[1]):
            if adj_mat[u, v] == 1 and growth_rate[v] > growth_rate[u]:
                num = growth_rate[v] - growth_rate[u]
                denom = 0
                for j in range(growth_rate.size):
                    if adj_mat[u, j] == 1 and growth_rate[j] > growth_rate[u]:
                        denom += (growth_rate[j] - growth_rate[u])
                subt_mat[u, v] = num / denom
        if subt_mat[u].sum() == 0:
            subt_mat[u, u] = 1
    # add immigration probabilities        
    if immigration_prob != 0:
        for i in range(subt_mat.shape[0]):
            num_non_zero = np.count_nonzero(subt_mat[i,:])
            num_zeros = subt_mat.shape[0] - num_non_zero
            for j in range(subt_mat.shape[1]):
                if subt_mat[i, j] == 0:
                    subt_mat[i, j] += immigration_prob / num_zeros
                else:
                    subt_mat[i, j] -= immigration_prob / num_non_zero

    return subt_mat

def seq_mat(ab_seq, immigration_prob = 0):
    M = np.identity(len(GENO))
    for ab in ab_seq:
        M_ab = cpm(ab, immigration_prob)
        M = np.dot(M, M_ab)
    return M


# change this to it takes in a list of transition matrices
def cycle_prob(transition_mat_list):
    '''
    Given a list of transition matrices corresponding to an antibiotic sequence
    and a probability model, returns the transition matrix for the antibiotic
    sequence. 
    '''
    M = np.identity(len(GENO))
    for mat in transition_mat_list:
        #files = os.listdir(ADJ_MAT_DIR)
        #ab_file = list(filter(lambda f: ab in re.findall('[A-Z]{2,3}', f), files))[0]
        #ab_path = '{}/{}/{}'.format(os.getcwd(), ADJ_MAT_DIR, ab_file)
        #M_ab = prob_model(ab_path)
        M = np.dot(M, mat)
    return M

def repeated_cycle_prob(seq_matrix, num_cycles):
    '''
    Caclulates the transition matrix for a given seqence of antibiotics after a given
    number of cycles 

    Inputs:
        seq_matrix: a numpy array of the transition matrix corresponding to an antibiotic sequence
        num_cycles: the number of cycles that the matrix will be repeated
        initial_state: if specified it should be the row vector containing the intial state probabilities
        plot: if True will generate a plot of the figure
        savefig: if not None, indicates the name of the file for the plot to be saved in
    '''
    M = np.identity(len(GENO))
    for i in range(num_cycles):
        M = np.dot(M, seq_matrix)
    return M

def eq_cycle_prob(mat):
    '''
    Finds the equilibrium (stationary) transition probabilities for a transition matrix
    '''
    M_eq = np.identity(len(GENO))
    count = 0
    while not np.all(np.dot(M_eq, mat) == M_eq) and count < 1000:
        M_eq = np.dot(M_eq, mat)
        count += 1
    return M_eq

def plot_cycle_probs(ab_seqs, initial_state, prob_model, legend = True, savefig = None, immigration_prob = 0):
    '''
    Creates a stacked bar plot of the distribution of genotypes after a specified
    number of cycles of an antibiotic sequence

    Inputs:
        seq_matrix: a numpy array of the transition matrix for a sequence of antibiotics
        initial_state: a numpy array of the initial probability distribution
        num_cycles: the number of cycles
        sasvefig: default value of None. If not none should be a string specifying the
            filename

    Returns:
        plot_info: a numpy array where the jth column contains the probability distribution
            after the jth cycle
    '''
    # plot info is an array containing probabilities where the rows are the genotype and 
    # the columns are the number of times the cycle has been repeated 

    num_cycles = min(50, len(ab_seqs))
    plot_info = np.zeros((len(GENO), num_cycles))
    M = np.identity(len(GENO))
    for i in range(num_cycles):
        #ab_seq = ab_seqs[i]
        transition_mat_list = []
        for ab in ab_seqs:
            transition_mat_list.append(prob_model(ab, immigration_prob))
        M = cycle_prob(transition_mat_list)
        M = np.dot(M, prob_model(ab_seqs[i], immigration_prob))
        initial_state = np.dot(initial_state, M)
        plot_info[:,i] = initial_state
    # create color palate     
    color = iter(plt.cm.cubehelix(np.linspace(0, 1, len(GENO)+3)))
    c = next(color)
    c = next(color)
    c = next(color)

    ind = np.arange(num_cycles)
    plt.bar(ind, plot_info[0,:num_cycles], color = c, label = '{}'.format(GENO[0]))
    old = copy.deepcopy(plot_info[0,:num_cycles])
    for row in range(1, plot_info.shape[0]):
        c = next(color)
        plt.bar(ind, plot_info[row,:num_cycles], bottom = old, color = c, label = '{}'.format(GENO[row]))
        old += plot_info[row,:num_cycles]
    plt.ylim(0,1)
    plt.ylabel('probability', fontsize = 9, fontname = "serif")

    plt.xlim(-0.2,num_cycles)
    plt.xlabel('antibiotic', fontsize = 9, fontname = "serif")
    
    #xlabels = ['{}'.format(num + 1) for num in range(num_cycles)]
    #xlabels = []
    #for ab_seq in ab_seqs:
    #    label = "-".join(list(ab_seq))
    #    xlabels.append(label)
    xlabels = ab_seqs
    if num_cycles > 20:
        plt.tick_params(axis='both', which='major', labelsize=9)
        plt.xticks(ind + 0.4, xlabels, rotation = 'vertical')
    else:
        plt.xticks(ind + 0.4, xlabels)
    
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    if legend:
        plt.legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, frameon = False, prop={"family":"serif", "size":9})
    
    if savefig != None:
        plt.savefig(savefig, bbox_inches = 'tight')
    #plt.close()
    print("final state?", plot_info[:,-1])
    return plot_info

def optimize_seq(k, initial_state, prob_model, ab_list = ANTIBIOTICS, eq = False, immigration_prob = 0):
    '''
    Finds the sequence of antibiotics of length k maximizing the probability of 
    returning to the susceptible genotype 
    '''
    ab_seqs = list(itertools.permutations(ab_list, k))
    max_prob = 0
    best_seq = []
    susc_index = GENO.index('0000')
    count = 1
    for seq in ab_seqs:
        if count % 100 == 0:
            print(count)
        count += 1
        seq_matrices = []
        for ab in seq:
            seq_matrices.append(prob_model(ab, immigration_prob))
        #prob = cycle_prob(seq, prob_model, initial_state = initial_state)[susc_index]
        M = cycle_prob(seq_matrices)
        if eq:
            M = eq_cycle_prob(M)
        prob = np.dot(initial_state, M)[susc_index]
        if max_prob < prob:
            max_prob = prob
            best_seq = [seq]
        elif max_prob == prob:
            best_seq.append(seq)
    return best_seq, max_prob

def adaptive_optimize(initial_state, prob_model, ab_list = ANTIBIOTICS, immigration_prob = 0, number = 0):
    #print("ANTIBIOTICS:", ANTIBIOTICS)
    growth_rate_cutoff = .001
    growth_rates = import_arr(GROWTH_RATE_FILE)
    #print("ab_list inside function:", ab_list)
    min_expected_growth_rate = float("inf")
    best_ab = None
    for ab in ab_list:
        M = prob_model(ab, immigration_prob)
        state = np.dot(initial_state, M)
        expected = 0
        ab_index = ANTIBIOTICS.index(ab)

        for i in range(len(state)):
            expected += growth_rates[ab_index, i] * state[i]
        if expected < min_expected_growth_rate:
            best_ab = ab
            min_expected_growth_rate = expected
    print("best ab", best_ab)
    print("expected gr", min_expected_growth_rate)
    return best_ab, min_expected_growth_rate

def re_optimize_seqs(k, initial_state, prob_model, num_cycles, immigration_prob = 0):
    '''
    Re-optimizes the sequences of length k after each completion of a sequence.
    The number of times this process is repeated is specified by num_cycles
    '''
    seqs = []
    best_seq, max_prob = optimize_seq(k, initial_state, prob_model, immigration_prob)
    #if there's a tie, pick the first sequence. 
    seqs.append(best_seq[0])
    state = cycle_prob(best_seq[0], prob_model, initial_state, immigration_prob)
    for i in range(1, num_cycles):
        print(i)
        best_seq, max_prob = optimize_seq(k, state, prob_model, immigration_prob)
        seqs.append(best_seq[0])
        state = cycle_prob(best_seq[0], prob_model, state, immigration_p)
    return seqs 

def generate_plots(ab_seqs, initial_state):
    '''
    Generates plots for ab_seqs given an initial state for 15 and 50 cycles
    '''
    for ab_seq in ab_seqs:
        repeated_cycle_prob(15, ab_seq, cpm, initial_state, plot = True, savefig = True)
        repeated_cycle_prob(50, ab_seq, cpm, initial_state, plot = True, savefig = True)
    return 

# given an antibiotic, find the average growth rate across all of the genotypes
def avg_growth_rate(ab):
    growth_rates = import_arr(GROWTH_RATE_FILE)

    ab_index = ANTIBIOTICS.index(ab)
    ab_rates = growth_rates[ab_index,:]
    avg = sum(ab_rates) / len(ab_rates)
    
    return avg

def adaptive_cycling(initial_state, prob_model, num_cycles, immigration_prob = 0, return_expected_gr = False, number = 0):
    
    growth_rate_cutoff = .001
    growth_rates = import_arr(GROWTH_RATE_FILE)
    #state = initial_state
    adaptive_seq = []
    ab_list = ['AMP','AM','CEC','CTX','ZOX','CXM','CRO','AMC', \
               'CAZ','CTT','SAM','CPR','CPD','TZP','FEP']
    ab_list.sort()
    prev_ab = None
    states = [initial_state]

    gr_list = []

    while len(adaptive_seq) < num_cycles:
        if prev_ab != None:
            ab_list.remove(prev_ab)
            ab_list.sort()
        #print("ab list", ab_list)
        best_ab, expected_gr = adaptive_optimize(initial_state, prob_model, ab_list = ab_list, immigration_prob = immigration_prob, number = number)
        #print("ab:", best_ab)
        
        gr_list.append(expected_gr)
        #best_ab, prob_susc = optimize_seq(1, state, prob_model, ab_list = ab_list, eq = False, immigration_prob = immigration_prob)
        #best_ab = best_ab[0][0]
        if prev_ab != None:
            ab_list.append(prev_ab)
            ab_list.sort()

        prev_ab = best_ab
        adaptive_seq.append(best_ab)
        #if number == 5:
            #print("best ab:", best_ab)
        M = prob_model(best_ab, immigration_prob)
        #M = np.dot(M, M_ab)
        initial_state = np.dot(initial_state, M)
        #print("state:", initial_state)
        states.append(initial_state)
        #print("expected_gr:", expected_gr)
        # find the expected growth rate
        ab_index = ANTIBIOTICS.index(best_ab)
        #expected = 0
        #print("state:", state)
        #for i in range(len(state)):
        #    expected += state[i] * growth_rates[ab_index, i]
       

        while expected_gr <= .001 and len(adaptive_seq) < num_cycles:
            adaptive_seq.append(best_ab)
            #print("ab:", best_ab)
            M = prob_model(best_ab, immigration_prob)
            initial_state = np.dot(initial_state, M)
            #print("state:", initial_state)
            states.append(initial_state)
            #print("state:", state)
            expected_gr = 0

            for i in range(len(initial_state)):
                #print("state[i]:", state[i])
                #print("growth_rates[ab_index, i]:", growth_rates[ab_index, i])
                #print("state[i] * growth_rates[ab_index, i]:", state[i] * growth_rates[ab_index, i])
                expected_gr += initial_state[i] * growth_rates[ab_index, i]
                #print("calculating expected_gr:", expected_gr)
            #print("EXPECTED:", expected_gr)
            gr_list.append(expected_gr)
            #print("expectetd_gr:", expected_gr)
    states_arr = np.array(states[-50:])
    #print(states_arr.shape)
    mean_state = np.mean(states_arr, axis = 0)
    #print("mean_state:", mean_state)
    #print(sum(mean_state))

    if return_expected_gr:
        return adaptive_seq, gr_list
    return adaptive_seq, mean_state

def random_cycling(initial_state, prob_model, num_cycles, immigration_prob = 0):
    growth_rate_cutoff = .001
    growth_rates = import_arr(GROWTH_RATE_FILE)
    state = initial_state
    adaptive_seq = []
    ab_list = ANTIBIOTICS
    states = [state]

    while len(adaptive_seq) < num_cycles:
        rand_int = random.randint(0, 14)
        rand_ab = ANTIBIOTICS[rand_int]
        adaptive_seq.append(rand_ab)
        #print("best ab:", best_ab)
        M = prob_model(rand_ab, immigration_prob)
        #M = np.dot(M, M_ab)
        state = np.dot(state, M)
        states.append(state)       
        # find the expected growth rate
        ab_index = ANTIBIOTICS.index(rand_ab)
        expected = 0
        for i in range(len(state)):
            expected += state[i] * growth_rates[ab_index, i]

        while expected <= growth_rate_cutoff and len(adaptive_seq) < num_cycles:
            adaptive_seq.append(rand_ab)
            state = np.dot(state, M)
            states.append(state)
            expected = 0
            for i in range(len(state)):
                expected += state[i] * growth_rates[ab_index, i]
    states_arr = np.array(states[-100:])
    #print(states_arr.shape)
    mean_state = np.mean(states_arr, axis = 0)
    #print(mean_state)
    #print(sum(mean_state))
    return adaptive_seq, mean_state


def plot_time_series(ab_seq, prob_model, num_cycles, initial_state, savefig = None, immigration_prob = 0):
    growth_rates = import_arr(GROWTH_RATE_FILE)
    growth_rate_cutoff = .001
    expected_rates = []
    count = 0
    xaxis = []
    xlabels = []
    #M = np.identity(16)
    for n in range(num_cycles):
        for ab in ab_seq:
            ab_index = ANTIBIOTICS.index(ab)
            M = prob_model(ab, immigration_prob)
            initial_state = np.dot(initial_state, M)
            val = 0
            for i in range(len(initial_state)):
                val += initial_state[i] * growth_rates[ab_index, i]
            expected_rates.append(val)
            xaxis.append(count)
            count += 1
            xlabels.append(ab)
    plt.plot(xaxis, expected_rates, "ro")
    plt.tick_params(axis='both', which='major', top = False, right = False, labelsize=8)
    plt.xticks(xaxis, xlabels, rotation = "vertical")
    plt.xlabel("antibiotic")
    plt.ylabel("expectetd growth rate")
    plt.xlim((-1, count+1))
    cycle = "-".join(ab_seq)

    if savefig != None:
        plt.savefig(savefig)
    plt.close()

    avg_growth_rate = sum(expected_rates[-100:]) / len(expected_rates[-100:])
    return avg_growth_rate

'''
if __name__ == '__main__':
    initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    
    for i in range(100):
        print(i)
        initial_state = []
        for j in range(16):
            initial_state.append(random.random())
        initial_state = np.array(initial_state)
        initial_state = initial_state / sum(initial_state)

        adaptive_seq = adaptive_cycling(initial_state, cpm, 50, immigration_prob = 0)
        plot_cycle_probs(adaptive_seq, initial_state, cpm, savefig = "random_initial_state/growth_rates/adaptive_seq_{}.png".format(i), immigration_prob = 0)
        adaptive_seq = adaptive_cycling(initial_state, cpm, 50, immigration_prob = 0.05)
        plot_cycle_probs(adaptive_seq, initial_state, cpm, savefig = "random_initial_state/growth_rates/adaptive_seq_imm_{}.png".format(i), immigration_prob = 0.05)
    
    adaptive_seq = adaptive_cycling(initial_state, cpm, 50, immigration_prob = 0)
    plot_cycle_probs(adaptive_seq, initial_state, cpm, savefig = "susc_prob_adaptive_seq.png", immigration_prob = 0)
    plot_time_series(adaptive_seq, cpm, 1, initial_state, savefig = "susc_prob_growth_rate_adaptive_time_series.png", immigration_prob = 0)
    adaptive_seq_imm = adaptive_cycling(initial_state, cpm, 50, immigration_prob = 0.05)
    plot_cycle_probs(adaptive_seq_imm, initial_state, cpm, savefig = "susc_prob_growth_rate_adaptive_seq_imm.png", immigration_prob = 0.05)
    plot_time_series(adaptive_seq_imm, cpm, 1, initial_state, savefig = "susc_prob_growth_rate_adaptive_time_series_imm.png", immigration_prob = 0.05)

    # optimal cycle determined by Mira et al.
    #plot_time_series(['AM', 'TZP'], cpm, 30, initial_state, immigration_prob = 0)
    #plot_time_series(['AM', 'TZP'], cpm, 30, initial_state, immigration_prob = 0.05)
    #plot_time_series(['AM', 'TZP'], cpm, 30, initial_state, expected = False, immigration_prob = 0)
    #plot_time_series(['AM', 'TZP'], cpm, 30, initial_state, expected = False, immigration_prob = 0.05)
    # optimal cycle considering equilibrium probabilities
    #plot_time_series(['CTX', 'CXM', 'TZP', 'CRO', 'AM'], cpm, 12, initial_state, immigration_prob = 0)
    #plot_time_series(['CTX', 'CXM', 'TZP', 'CRO', 'AM'], cpm, 12, initial_state, immigration_prob = 0.05)
    #plot_time_series(['CTX', 'CXM', 'TZP', 'CRO', 'AM'], cpm, 12, initial_state, expected = False, immigration_prob = 0)
    #plot_time_series(['CTX', 'CXM', 'TZP', 'CRO', 'AM'], cpm, 12, initial_state, expected = False, immigration_prob = 0.05)
    ab_seq = ["AM", "TZP"]

    M = seq_mat(ab_seq)
    ## Short Period ##
    
    #Single cycle
    p = np.dot(initial_state, M)[0]
    print("Short period single Cycle p:", p)
    #Infinite cycle
    p = np.dot(initial_state, eq_cycle_prob(M))[0]
    print("Short period Infinite Cycle p:", p)

    ## Long period ##

    M = np.identity(len(GENO))
    for ab in ab_seq:
        M_ab = eq_cycle_prob(cpm(ab))
        M = np.dot(M, M_ab)

    #Single cycle
    p = np.dot(initial_state, M)[0]
    print("Long period-single cycle p", p)
    #Infinite cycle
    p = np.dot(initial_state, eq_cycle_prob(M))[0]
    print("Long period -infinie cycle p:", p)
    
    # prepare growth rates table
    f = open("growth_rate.csv", "r")
    f_new = open("growth_rate_table.csv", "w")
    geno_str = [str(geno) for geno in GENO]
    f_new.write("," + ",".join(geno_str) + "\n")
    ab_index = 0
    for line in f:
        line_list = line.strip().split(",")
        for i in range(len(line_list)):
            line_list[i] = str(round(float(line_list[i]) * 1000, 3))
        new_line = ",".join(line_list)
        f_new.write(ANTIBIOTICS[ab_index] + "," + new_line + "\n") 
        ab_index += 1
    f.close()
    f_new.close()
    '''
    

