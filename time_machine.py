import itertools
import random
import numpy as np
import pandas as pd
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
    while not np.all(np.dot(M_eq, mat) == M_eq) and count < 100000:
        M_eq = np.dot(M_eq, mat)
        count += 1
    return M_eq

def plot_cycle_probs(seq_matrix, initial_state, num_cycles, savefig = None):
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
    plot_info = np.zeros((seq_matrix.shape[0], num_cycles))
    plot_info[:,0] = np.dot(initial_state, seq_matrix)
    M = np.identity(len(GENO))
    for i in range(1, num_cycles):
        M = np.dot(M, seq_matrix)
        plot_info[:,i] = np.dot(initial_state, M)

    # create color palate     
    color = iter(plt.cm.cubehelix(np.linspace(0, 1, len(GENO)+3)))
    c = next(color)
    c = next(color)
    c = next(color)

    ind = np.arange(num_cycles)
    plt.bar(ind, plot_info[0,:], color = c, label = '{}'.format(GENO[0]))
    old = plot_info[0,:]
    for row in range(1, plot_info.shape[0]):
        c = next(color)
        plt.bar(ind, plot_info[row,:], bottom = old, color = c, label = '{}'.format(GENO[row]))
        old += plot_info[row,:]
    plt.ylim(0,1)
    plt.ylabel('Probability')
    plt.xlim(-0.2,num_cycles)
    plt.xlabel('Cycle')
    xlabels = ['{}'.format(num + 1) for num in range(num_cycles)]
    
    if num_cycles > 20:
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.xticks(ind + 0.4, xlabels, rotation = 'vertical')
    else:
        plt.xticks(ind + 0.4, xlabels)
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, fontsize='small')
    
    if savefig != None:
        plt.savefig(savefig, bbox_inches = 'tight')
        plt.close()
    else:
        plt.show()

    return plot_info 

def generate_latex_tables(ab_seqs, initial_state):
    '''
    This function just formats the output probabilities for repeated_cycle_prob into 
    a rows of a latex tabel
    '''
    for ab_seq in ab_seqs:
        file_name = '_'.join(ab_seq)
        f = open('tables/{}'.format(file_name), 'w')
        for cycle_num in range(1, 16):
            M_row = repeated_cycle_prob(cycle_num, ab_seq, cpm, initial_state)
            probs = '&'.join([str(x) for x in M_row])
            f.write('Cycle {},{} {}\n'.format(cycle_num, probs, '\\\\'))
    f.close()
    return

def optimize_seq(k, initial_state, prob_model, eq = False):
    '''
    Finds the sequence of antibiotics of length k optimizing the probability of 
    returning to the susceptible genotype 
    '''
    ab_seqs = list(itertools.permutations(ANTIBIOTICS, k))
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
            seq_matrices.append(prob_model(ab))
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

def re_optimize_seqs(k, initial_state, prob_model, num_cycles):
    '''
    Re-optimizes the sequences of length k after each completion of a sequence.
    The number of times this process is repeated is specified by num_cycles
    '''
    seqs = []
    best_seq, max_prob = optimize_seq(k, initial_state, prob_model)
    #if there's a tie, pick the first sequence. 
    seqs.append(best_seq[0])
    state = cycle_prob(best_seq[0], prob_model, initial_state)
    for i in range(1, num_cycles):
        print(i)
        best_seq, max_prob = optimize_seq(k, state, prob_model)
        seqs.append(best_seq[0])
        state = cycle_prob(best_seq[0], prob_model, state)
    return seqs 

def generate_plots(ab_seqs, initial_state):
    '''
    Generates plots for ab_seqs given an initial state for 15 and 50 cycles
    '''
    for ab_seq in ab_seqs:
        repeated_cycle_prob(15, ab_seq, cpm, initial_state, plot = True, savefig = True)
        repeated_cycle_prob(50, ab_seq, cpm, initial_state, plot = True, savefig = True)
    return 

def compare(ab_seq, prob_model, initial_states, num_cycles, Mira = True, Nichol = True, immigration_prob = 0):
    '''
    This function should be run in a directory named the sequence of antibiotics
    '''
    if Mira:
        mat_list_Mira = []
    if Nichol:
        mat_list_Nichol = []
    for ab in ab_seq:
        M_ab = prob_model(ab, immigration_prob = immigration_prob)
        if Mira:
            mat_list_Mira.append(M_ab)
        if Nichol:
            mat_list_Nichol.append(eq_cycle_prob(M_ab))

    if Mira:
        seq_mat_Mira = cycle_prob(mat_list_Mira)
        for state in initial_states:
            if np.all(state == np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])):
                initial = 'initial_state_susceptible'
            elif np.all(state == np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])):
                initial = 'initial_state_resistant'
            elif np.all(state == np.array([1/16] * 16)):
                initial = 'initial_state_unknown'
            else:
                initial = '{}'.format(state)
            filename = 'new_plots/CEC_SAM_CTX_AM_{}_{}_{}_immigration={}.png'.format(num_cycles, prob_model.__name__, initial, immigration_prob)
            #plot_cycle_probs(seq_mat_Mira, state, num_cycles, savefig = filename)
            plot_cycle_probs(seq_mat_Mira, state, num_cycles, savefig = filename)
    if Nichol:
        seq_mat_Nichol = cycle_prob(mat_list_Nichol)
        for state in initial_states:
            if np.all(state == np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])):
                initial = 'initial_state_susceptible'
            elif np.all(state == np.array([1/16] * 16)):
                initial = 'initial_state_unknown'
            else:
                initial = '{}'.format(state)
            filename_15 = 'new_plots/Nichol_{}_{}_{}_immigration={}.png'.format(num_cycles, prob_model.__name__, initial, immigration_prob)
            plot_cycle_probs(seq_mat_Nichol, state, num_cycles, savefig = filename_15)

    return

if __name__ == '__main__':
    initial_state = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    #compare(['AM', 'TZP'], cpm, [initial_state], 15, Nichol = False)
    #compare(['AM', 'TZP'], cpm, [initial_state], 15, Nichol = False, immigration_prob = 0.05)
    #compare(['CEC', 'CTX', 'TZP', 'CTX', 'TZP', 'AM'], cpm, [initial_state], 15, Nichol = False)
    compare(['CEC', 'SAM', 'CTX', 'AM'], cpm, [initial_state], 15, Nichol = False)
