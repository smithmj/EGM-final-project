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

ADJ_MAT_DIR = 'adj_mat'

GROWTH_RATE_FILE = 'growth_rate.csv'

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

# Gives the average growth rates. Column labels are given by geno array.
# Rows labels are given by antibiotics array. Each row is a fitness landscape
# of the genotype under the given antibiotic. 
avg_growth = import_arr('growth_rate.csv')

def epm(adj_mat_file):
    '''
    Creates the transition matrix for a given adjacency matrix 
    based on the equal probability model
    '''
    subt_mat = import_arr(adj_mat_file)
    for i in range(subt_mat.shape[0]):
        s = subt_mat[i].sum()
        if s != 0:
            subt_mat[i] = subt_mat[i] / subt_mat[i].sum()
    return subt_mat

def cpm(adj_mat_file, growth_rate_file):
    '''
    Creates the transition matrix for a given adjacency matrix 
    based on the correlated probability model
    '''
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
    return subt_mat

def cycle_prob(ab_seq, prob_model, initial_state = None):
    '''
    Given an antibiotic sequence, the index of the initial genotype, a probability
    model, adjacency matrix filename, and the growth rate filename, returns the probability
    of ending with a susceptible genotype. 
    '''
    M = np.identity(16)
    for ab in ab_seq:
        files = os.listdir(ADJ_MAT_DIR)
        ab_file = list(filter(lambda f: ab in re.findall('[A-Z]{2,3}', f), files))[0]
        ab_path = '{}/{}/{}'.format(os.getcwd(), ADJ_MAT_DIR, ab_file)
        M_ab = prob_model(ab_path, GROWTH_RATE_FILE)
        M = np.dot(M, M_ab)
    if initial_state != None:
        return np.dot(initial_state, M)
    return M

def repeated_cycle_prob(ab_seqs, prob_model, initial_state = None, plot = False, savefig = None):
    '''
    Caclulates and potentially plots the probability for a given seqence of antibiotics 

    Inputs:
        ab_seqs: a list of tuples of antibiotic sequences. For example: [('AM', 'TZP')] * 15 would give 
            fifteen cycles of ('AM', 'TZP') sequence.
        prob_model: either cpm or epm
        initial_state: if specified it should be the row vector containing the intial state probabilities
        plot: if True will generate a plot of the figure
        savefig: if not None, indicates the name of the file for the plot to be saved in


    '''
    num_cycles = len(ab_seqs)
    M = cycle_prob(ab_seqs[0], prob_model)
    if plot:
        plot_info = np.zeros((M.shape[0], num_cycles))
        plot_info[:,0] = np.dot(initial_state, M)
    for i in range(1, num_cycles):
        M = np.dot(M, cycle_prob(ab_seqs[i], prob_model))
        if plot:
            plot_info[:,i] = np.dot(initial_state, M)

    if plot:
        all_zeros = np.zeros((num_cycles))
        ind = np.arange(num_cycles)
        n = 0
        for i in range(plot_info.shape[0]):
            if not np.all(plot_info[i,:] == all_zeros):
                n += 1
        
        color = iter(plt.cm.cubehelix(np.linspace(0,1,n+2)))
        c = next(color)
        c = next(color)
        row = 0
        while np.all(plot_info[row,:] == all_zeros):
            row += 1
        plt.bar(ind, plot_info[row,:], color = c, label = '{}'.format(GENO[0]))
        old = plot_info[row,:]
        for row in range(1, plot_info.shape[0]):
            if not np.all(plot_info[row,:] == all_zeros):
                c = next(color)
                plt.bar(ind, plot_info[row,:], bottom = old, color = c, label = '{}'.format(GENO[row]))
                old += plot_info[row,:]
        plt.ylim(0,1)
        plt.ylabel('Probability')
        plt.xlim(-0.2,num_cycles)
        plt.xlabel('Cycle')
        #xlabels = ['{}'.format(num + 1) for num in range(num_cycles)]
        xlabels = ['{}'.format('-'.join(seq)) for seq in ab_seqs]
        if num_cycles > 20:
            plt.tick_params(axis='both', which='major', labelsize=10)
            plt.xticks(ind + 0.4, xlabels, rotation = 'vertical')
        else:
            plt.xticks(ind + 0.4, xlabels)
        plt.margins(0.2)
        plt.subplots_adjust(bottom=0.15)
        plt.legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, fontsize='small')
        if savefig != None:
            f = '_'.join(ab_seq)
            plt.savefig(savefig, bbox_inches = 'tight')
            plt.close()
        else:
            plt.show()

    if initial_state != None:
        return np.dot(initial_state, M)
    return M

def generate_plots(ab_seqs, initial_state):
    '''
    Generates plots for ab_seqs given an initial state for 15 and 50 cycles
    '''
    for ab_seq in ab_seqs:
        repeated_cycle_prob(15, ab_seq, cpm, initial_state, plot = True, savefig = True)
        repeated_cycle_prob(50, ab_seq, cpm, initial_state, plot = True, savefig = True)
    return 

def generate_tables(ab_seqs, initial_state):
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
            f.write('Cycle {},{} {}\n'.format(cycle_num, probs, '\\'))
    f.close()
    return

def optimize_seq(k, initial_state, prob_model):
    '''
    Finds the sequence of antibiotics of length k optimizing the probability of 
    returning to the susceptible genotype 
    '''
    #ab_seqs = list(itertools.product(ANTIBIOTICS, repeat = k))
    ab_seqs = list(itertools.permutations(ANTIBIOTICS, k))
    max_prob = 0
    best_seq = []
    susc_index = GENO.index('0000')
    for seq in ab_seqs:
        prob = cycle_prob(seq, prob_model, initial_state = initial_state)[susc_index]
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












