import numpy as np
import matplotlib.pyplot as plt
import heapq
import random
import copy
import pandas as pd
import scipy.stats
import time_machine
import rankings
import correlation

if __name__ == "__main__":
    
    # S1 (Example plots for long term cycling effects)
    '''
    initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    fig, axarr = plt.subplots(1, 2)
    n = -1
    for ab_seq in [['AM', 'TZP'], ['CEC', 'SAM', 'CTX', 'AM']]:
        n += 1
        num_cycles = 15
        plot_info = np.zeros((len(time_machine.GENO), num_cycles))
        transition_mat_list = []
        for ab in ab_seq:
            transition_mat_list.append(time_machine.cpm(ab, immigration_prob = 0))
        M_seq = time_machine.cycle_prob(transition_mat_list)
        M = np.identity(len(time_machine.GENO))
        for j in range(num_cycles):
            #M = np.dot(M_seq, time_machine.cpm(ab_seq[j], immigration_prob = 0))
            M = np.dot(M, M_seq)
            initial_state = np.dot(initial_state, M)
            plot_info[:,j] = initial_state     
        color = iter(plt.cm.cubehelix(np.linspace(0, 1, len(time_machine.GENO)+3)))
        c = next(color)
        c = next(color)
        c = next(color)

        ind = np.arange(num_cycles)
        axarr[n].bar(ind, plot_info[0,:num_cycles], color = c, label = '{}'.format(time_machine.GENO[0]))
        old = copy.deepcopy(plot_info[0,:num_cycles])
        for row in range(1, plot_info.shape[0]):
            c = next(color)
            axarr[n].bar(ind, plot_info[row,:num_cycles], bottom = old, color = c, label = '{}'.format(time_machine.GENO[row]))
            old += plot_info[row,:num_cycles]
        axarr[n].set_ylim(0,1)
        axarr[n].tick_params(axis='both', which='major', labelsize=9)
        if n == 0:
            axarr[n].set_ylabel('frequency', fontsize = 9, fontname = "serif")
            #axarr[n].tick_params(fontsize = 9, fontname = "serif")


        axarr[n].set_xlim(-0.2,num_cycles)
        axarr[n].set_xlabel("{} cycle".format('-'.join(ab_seq)), fontsize = 9, fontname = "serif")
        #axarr[n].set_xticklabels(ind, ind /2+ 0.2, fontsize = 9, fontname = "serif")
        #axarr[n].get_xaxis().set_ticks([])
        if n == 1:
            axarr[n].legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, frameon = False, prop={"family":"serif", "size":9})
    
            #axarr[n].margins(0.2)
            #axarr[n].subplots_adjust(bottom=0.15)
        
        
    #fig.tight_layout()
    #fig.subplots_adjust(top = 0.1)
    fig.set_size_inches(7.28, 4, forward = True)
    
    plt.savefig("figures/supplement/repeated_cycle_probs.pdf")
    plt.close()
    '''
    
    #S4 (Adaptive cycling starting from 0000)
    
    fig, axarr = plt.subplots(2 ,2)
    
    #plot_cycle_probs(adaptive_seq, initial_state, cpm, savefig = "susc_prob_adaptive_seq.png", immigration_prob = 0)
    #for ip in [0, 0.01]:
    for ip in [0, .01]:
        print("ip:", ip)
        if ip == 0:
            n = 0
        elif ip == 0.01:
            n = 1
        initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        adaptive_seq, expected_gr = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 30, immigration_prob = ip, return_expected_gr = True)
        print("adaptive_seq", adaptive_seq)
        #print("expected_gr", expected_gr)
        # plot cycle probs
        num_cycles = len(adaptive_seq)
        plot_info = np.zeros((len(time_machine.GENO), num_cycles))
        M = np.identity(len(time_machine.GENO))
        
        for j in range(len(adaptive_seq)):
            ab = adaptive_seq[j]
            M = np.dot(M, time_machine.cpm(ab, immigration_prob = ip))
            initial_state = np.dot(initial_state, M)
            plot_info[:,j] = initial_state 

        color = iter(plt.cm.cubehelix(np.linspace(0, 1, len(time_machine.GENO)+3)))
        c = next(color)
        c = next(color)
        c = next(color)

        ind = np.arange(num_cycles)
        axarr[1,n].bar(ind, plot_info[0,:], color = c, label = '{}'.format(time_machine.GENO[0]))
        old = copy.deepcopy(plot_info[0,:])
        for row in range(1, plot_info.shape[0]):
            c = next(color)
            axarr[1,n].bar(ind, plot_info[row,:], bottom = old, color = c, label = '{}'.format(time_machine.GENO[row]))
            old += plot_info[row,:]
        axarr[1,n].set_ylim(0,1)
        if n == 0:
            axarr[1,n].set_ylabel('frequency', fontsize = 9, fontname = "serif")
        axarr[1,n].set_xlim(-0.2,num_cycles)
        axarr[1,n].set_xlabel('step', fontsize = 9, fontname = "serif")
        axarr[1,n].tick_params(axis='both', which='major', labelsize=9)
        axarr[1,n].set_xticks(ind + 0.4)
        axarr[1, n].set_xticklabels(adaptive_seq, fontsize = 6, fontname = "serif", rotation = 'vertical')
        #axarr[n].get_xaxis().set_ticks([])
        #axarr[n].set_yticklabels(fontsize = 9, fontname = "serif")
        axarr[1,n].margins(0.2)
        #axarr[n].subplots_adjust(bottom=0.15)
        if n == 1:
            axarr[1,n].legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, frameon = False, prop={"family":"serif", "size":9})

        #plot time series 
        #plot_time_series(adaptive_seq, cpm, 1, initial_state, savefig = "susc_prob_growth_rate_adaptive_time_series.png", immigration_prob = 0)
        growth_rates = time_machine.import_arr(time_machine.GROWTH_RATE_FILE)
        growth_rate_cutoff = .001
        expected_rates = []
        count = 0
        xaxis = []
        xlabels = []
        
        initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        for j in range(len(adaptive_seq)):
            ab = adaptive_seq[j]
            ab_index = time_machine.ANTIBIOTICS.index(ab)
            M = time_machine.cpm(ab, immigration_prob = ip)
            initial_state = np.dot(initial_state, M)
            val = 0
            for k in range(len(initial_state)):
                val += initial_state[k] * growth_rates[ab_index, k]
                
            expected_rates.append(val)
            xaxis.append(count)
            count += 1
        #print("expected_rates:", expected_rates)
        #print(len(expected_rates))
        #axarr[0,n].plot(xaxis, expected_rates, "ro", color = "gray")
        axarr[0,n].plot(xaxis, expected_gr, "ro", color = "gray")
        axarr[0,n].tick_params(axis='both', which='major', top = False, right = False, labelsize=9)
        
        axarr[0,n].get_xaxis().set_ticks([])
        #axarr[1,n].set_xticklabels(xlabels, xlabels, rotation = "vertical")
        
        #axarr[0,n].set_xlabel("step", fontsize = 9, fontname = "serif")
        if n == 0:
            axarr[0,n].set_ylabel("expectetd growth rate", fontsize = 9, fontname = "serif")
        axarr[0,n].set_xlim((-1, count+1))

    #fig.tight_layout()
    #fig.subplots_adjust(left = 0.2)
    fig.set_size_inches(7.28, 8, forward = True)
    
    plt.savefig("figures/supplement/adaptive_cycling_from_0000_new.pdf")
    plt.close()



    #adaptive_seq_imm = time_machine.adaptive_cycling(initial_state, cpm, 50, immigration_prob = 0.01)
    #plot_cycle_probs(adaptive_seq_imm, initial_state, cpm, savefig = "susc_prob_growth_rate_adaptive_seq_imm.png", immigration_prob = 0.01)
    #plot_time_series(adaptive_seq_imm, cpm, 1, initial_state, savefig = "susc_prob_growth_rate_adaptive_time_series_imm.png", immigration_prob = 0.01)

    '''
    #S5
    intial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    cycle_dict = {}
    g = open("figures/figure_4/random_initial_states_imm.txt", "r")
    #for i in range(100):
    i = 0
    fig, axarr = plt.subplots(1, 3)
    xlabels = [0, 5, 10, 15, 20, 25, 30]
    plt.setp(axarr, xticks=[0, 1/6*30, 2/6*30, 3/6*30, 4/6*30, 5/6*30, 30], xticklabels=xlabels)
    for line in g:
        i += 1
        print(i)
        l = line.strip().split(",")
        initial_state = np.array(l).astype("float")
    
        adaptive_seq = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 150, immigration_prob = 0.01)
        pattern = pattern_finder(adaptive_seq, min_reps = 10)
        len_pattern = len(pattern)
        pattern = "-".join(pattern)
        if i == 18 or i == 21 or i == 81:
            if i == 18:
                n = 0
                legend = False
            elif i == 21:
                n = 1
                legend = False
            elif i == 81:
                n = 2
                legend = True
            num_cycles = 30
            print(adaptive_seq)
            plot_info = np.zeros((len(time_machine.GENO), num_cycles))
            M = np.identity(len(time_machine.GENO))
            for j in range(num_cycles):
                transition_mat_list = []
                for ab in adaptive_seq:
                    transition_mat_list.append(time_machine.cpm(ab, immigration_prob = 0.01))
                M = time_machine.cycle_prob(transition_mat_list)
                M = np.dot(M, time_machine.cpm(adaptive_seq[j], immigration_prob = 0.01))
                initial_state = np.dot(initial_state, M)
                plot_info[:,j] = initial_state     
            color = iter(plt.cm.cubehelix(np.linspace(0, 1, len(time_machine.GENO)+3)))
            c = next(color)
            c = next(color)
            c = next(color)

            ind = np.arange(num_cycles)
            axarr[n].bar(ind, plot_info[0,:num_cycles], color = c, label = '{}'.format(time_machine.GENO[0]))
            old = copy.deepcopy(plot_info[0,:num_cycles])
            for row in range(1, plot_info.shape[0]):
                c = next(color)
                axarr[n].bar(ind, plot_info[row,:num_cycles], bottom = old, color = c, label = '{}'.format(time_machine.GENO[row]))
                old += plot_info[row,:num_cycles]
            axarr[n].set_ylim(0,1)
            if n == 0:
                axarr[n].set_ylabel('frequency', fontsize = 9, fontname = "serif")

            axarr[n].set_xlim(-0.2,num_cycles)
            if n == 1:
                axarr[n].set_xlabel('step', fontsize = 9, fontname = "serif")

            xlabels = adaptive_seq
            print(xlabels)
            axarr[n].tick_params(axis='both', which='major', labelsize=9)
            #axarr[n].set_xticklabels(xlabels, ind + 0.4, fontsize = 9, fontname = "serif", rotation = 'vertical')
            #axarr[n].get_xaxis().set_ticks([])
            #axarr[n].set_yticklabels(fontsize = 9, fontname = "serif")
    
            axarr[n].margins(0.2)
            #axarr[n].subplots_adjust(bottom=0.15)
            if legend:
                axarr[n].legend(bbox_to_anchor=(1., 1),loc=2, borderaxespad=0, frameon = False, prop={"family":"serif", "size":9})
            if n == 2:
                break
    #fig.tight_layout()
    fig.subplots_adjust(left = 0.06)
    fig.set_size_inches(7.28, 4, forward = True)
    
    plt.savefig("figures/supplement/no_pattern_plot.pdf")
    plt.close()
    '''
    

    #Supplementary table 2
    '''
    f = open("growth_rate_table.csv", "r")
    g = open("new_growth_rate_table.csv", "w")
    header = f.readline()
    g.write(header + ",mean\n")
    for line in f:
        l = line.split(',')
    '''