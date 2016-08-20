import numpy as np
#import matplotlib as mpl 
#mpl.rc('font',family='serif')
import matplotlib.pyplot as plt
import heapq
import random
import copy
import pandas as pd
import scipy.stats
import time_machine
import rankings
import correlation

CYCLES = ['FEP', 'AMC', 'SAM', 'FEP-SAM', 'FEP-CPD-AMC-SAM', 'SAM-AMC-FEP-CPD', 'FEP-CPD', 'CPD-SAM', 'AMC-SAM-CPD-FEP', \
          'CAZ', 'FEP-AMC-CPD-SAM', 'CPR', 'AMC-CPD-FEP-SAM', 'CPD-AMC-FEP-SAM', 'CTT', 'CAZ-CPR-CPD-AMC-TZP-SAM-CTT-FEP', \
          'AMC-CTT-SAM-CPR-CAZ-CPD-FEP-TZP', 'CPR-CAZ-AMC-TZP-SAM-FEP-CTT-CPD', 'CPR-CPD-AMC-TZP-SAM-FEP-CTT-CAZ', \
          'CTT-TZP-AMC-FEP-SAM-CPD-CAZ-CPR', 'AMC-FEP-CPR-CAZ-CTT-SAM-TZP-CPD', 'CAZ-CPR-CPD-CTT-TZP-SAM-FEP-AMC', \
          'CAZ-SAM-TZP-AMC-FEP-CTT-CPD-CPR', 'CPD', 'CPD-AMC-FEP-CAZ-CPD-AMC-FEP-CPR', 'CPD-AMC-FEP-CAZ-CPD-AMC-FEP-CPR', \
          'CPD-AMC-SAM-CAZ-CPD-AMC-SAM-CPR', 'CPD-CPR-AMC-TZP-SAM-FEP-CTT-CAZ', 'CPD-SAM-FEP-CPR-CPD-SAM-FEP-CAZ', \
          'CPR-CAZ-FEP-TZP-AMC-SAM-CTT-CPD', 'CPR-CAZ-SAM-TZP-AMC-CTT-CPD-FEP', 'CPR-CAZ-TZP-AMC-FEP-SAM-CTT-CPD', \
          'CPR-FEP-AMC-CAZ-TZP-CPD-SAM-CTT', 'CPR-SAM-AMC-CPD-CAZ-SAM-AMC-CPD', 'CTT-TZP-AMC-SAM-CPD-FEP-CAZ-CPR', \
          'FEP-CPD-CAZ-CPR-SAM-TZP-AMC-CTT', 'FEP-CPR-CAZ-CTT-TZP-AMC-SAM-CPD', 'FEP-CPR-TZP-AMC-CTT-SAM-CAZ-CPD', \
          'FEP-TZP-CAZ', 'SAM-CTT-CPR-TZP-CPD-CAZ-AMC-FEP', 'TZP-SAM-FEP-AMC-CPD-CAZ-CPR-CTT']

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        print("height", height)
        ax.text(rect.get_x() + rect.get_width()/2., height+.005,
                '{}'.format(round(height,2)),ha='center', va='bottom', fontname = "serif", fontsize = 9)
    return

def pattern_finder(ab_list, min_reps = 50):
    d = {}
    for l in range(1, 21):
        if l not in d:
            d[l] = ab_list[-l:]
        for rep in range(2, min_reps + 1):
            seq = ab_list[-l*rep:-l*(rep-1)]
            if seq != d[l]:
                d[l] = []
                break
    pattern_lens = [l for l in d if d[l] != []]
    if pattern_lens == []:
        return ["None"]
    pattern_len = min(pattern_lens)
    return d[pattern_len]

def is_pattern_eq(pattern1, pattern2):
    eq = False
    if pattern1 == pattern2:
        eq = True
    else:
        for i in range(len(pattern1)):
            first = pattern1.pop(0)
            pattern1.append(first)
            if pattern1 == pattern2:
                eq = True
                break
    return eq


if __name__ == "__main__":
    
    # make files
    '''
    for cycle in CYCLES:
        f = open("figures/figure_3/cycles/{}.txt".format(cycle), "w")
        f.close()
    '''
    ## Figure 1 ##
    '''
    f, axarr = plt.subplots(3, 1)

    rank1 = []
    rank2 = []
    f = open("figures/figure_1/single_cycle_vs_inf_cycle_ranks.txt", "r")
    for line in f:
        line_list = line.split(",")
        rank1.append(int(line_list[0].strip()))
        rank2.append(int(line_list[1].strip()))

    axarr[0].plot(rank1, rank2, 'ro', color = "gray")
    axarr[0].set_xlabel('single cycle\n(short period) rank', fontsize = 9, fontname = "serif")
    axarr[0].set_ylabel('infinite cycle\n(short period) rank', fontsize = 9, fontname = "serif")
    axarr[0].tick_params(axis='both', which='major', labelsize=9)
    axarr[0].tick_params(axis='both', which='minor', labelsize=9)
    max_x = max(rank1) + 2
    max_y = max(rank2) + 1000
    axarr[0].set_xlim(-1, max_x)
    axarr[0].set_ylim(-1000, max_y)
    f.close()

    rank3 = []
    rank4 = []
    f = open("figures/figure_1/inf_cycle_vs_long_period_ranks.txt", "r")
    for line in f:
        line_list = line.split(",")
        rank3.append(int(line_list[0].strip()))
        rank4.append(int(line_list[1].strip()))
        
    axarr[1].plot(rank3, rank4, 'ro', color = "gray")
    axarr[1].set_xlabel('infinite cycle\n(short period) rank', fontsize = 9, fontname = "serif")
    axarr[1].set_ylabel('single cycle\n(long period) rank', fontsize = 9, fontname = "serif")
    axarr[1].tick_params(axis='both', which='major', labelsize=9)
    axarr[1].tick_params(axis='both', which='minor', labelsize=9)
    max_x = max(rank3) + 2
    max_y = max(rank4) + 1000
    axarr[1].set_xlim(-1, max_x)
    axarr[1].set_ylim(-1000, max_y)
    f.close()

    rank5 = []
    rank6 = []
    f = open("figures/figure_1/short_period_vs_long_period_ranks.txt", "r")
    for line in f:
        line_list = line.split(",")
        rank5.append(int(line_list[0].strip()))
        rank6.append(int(line_list[1].strip()))
    f.close()
    
    axarr[2].plot(rank5, rank6, 'ro', color = "gray")
    axarr[2].set_xlabel('short period\n(single cycle) rank', fontsize = 9, fontname = "serif")
    axarr[2].set_ylabel('long period\n(single cycle) rank', fontsize = 9, fontname = "serif")
    axarr[2].tick_params(axis='both', which='major', labelsize=9)
    axarr[2].tick_params(axis='both', which='minor', labelsize=9)
    max_x = max(rank5) + 2
    max_y = max(rank6) + 1000
    axarr[2].set_xlim(-1, max_x)
    axarr[2].set_ylim(-1000, max_y)

    fig = plt.gcf()
    fig.tight_layout()
    fig.subplots_adjust(left=0.25)
    fig.set_size_inches(3.54, 7, forward = True)
    
    plt.savefig("figures/figure_1/rank_plots_vertical.pdf")
    
    '''    
    '''
    ## Figure 2 ##
    initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    #single_no_imm = []
    #single_with_imm = []
    #inf_no_imm = []
    #inf_with_imm = []
    with_imm = []
    no_imm = []
    
    #seqs = []
    labels = []

    print("1")
    print("s-c, s-p")
    (prob1, seq1) = rankings.combine_rankings("rankings/results/seq_ranks_2.txt", "rankings/results/seq_ranks_3.txt", \
                                                "rankings/results/seq_ranks_4.txt", "rankings/results/seq_ranks_5.txt", \
                                                "rankings/results/seq_ranks.txt")
    print(prob1, seq1)
    seq1 = seq1.split(",")
    seq1 = tuple([ab.strip("(").strip().strip(")").strip().strip("'") for ab in seq1])
    # single cycle without immigration
    mat1_no_imm = time_machine.seq_mat(seq1)
    p_single_no_imm = np.dot(initial_state, mat1_no_imm)[0]
    print(p_single_no_imm)
    #no_imm.append(p_single_no_imm)
    no_imm.append(prob1)    
    #single_no_imm.append(p_single_no_imm)
    # single cycle with immigration
    mat1_with_imm = time_machine.seq_mat(seq1, immigration_prob = 0.01)
    p_single_with_imm = np.dot(initial_state, mat1_with_imm)[0]
    with_imm.append(p_single_with_imm)
    #single_with_imm.append(p_single_with_imm)
    # infinite cycles without immigration
    #mat1_inf_no_imm = time_machine.eq_cycle_prob(mat1_no_imm)
    #p_inf_no_imm = np.dot(initial_state, mat1_inf_no_imm)[0]
    #inf_no_imm.append(p_inf_no_imm)
    #no_imm.append(p_inf_no_imm)
    # infinite cycles with immigration
    #mat1_inf_with_imm = time_machine.eq_cycle_prob(mat1_with_imm)
    #p_inf_with_imm = np.dot(initial_state, mat1_inf_with_imm)[0]
    #inf_with_imm.append(p_inf_with_imm)
    #with_imm.append(p_inf_with_imm)
    # sequence
    #seqs.append("-".join(seq1))
    labels.append("s-c\ns-p")
    

    print("2")
    print("s-c, l-p")
    (prob3, seq3) = rankings.combine_rankings("rankings/results/long_period_ranks_2.txt", "rankings/results/long_period_ranks_3.txt", \
                                                 "rankings/results/long_period_ranks_4.txt", "rankings/results/long_period_ranks_5.txt", \
                                                 "rankings/results/long_period_ranks.txt")
    print(prob3, seq3)
    seq3 = seq3.split(",")
    seq3 = tuple([ab.strip("(").strip().strip(")").strip().strip("'") for ab in seq3])
    # single cycle without immigration
    mat3_no_imm = time_machine.seq_mat(seq3)
    p_single_no_imm = np.dot(initial_state, mat3_no_imm)[0]
    print(p_single_no_imm)
    #no_imm.append(p_single_no_imm)
    no_imm.append(prob3)
    # single cycle with immigration
    mat3_with_imm = time_machine.seq_mat(seq3, immigration_prob = 0.01)
    p_single_with_imm = np.dot(initial_state, mat3_with_imm)[0]
    with_imm.append(p_single_with_imm)
    # infinite cycles without immigration
    #mat3_inf_no_imm = time_machine.eq_cycle_prob(mat3_no_imm)
    #p_inf_no_imm = np.dot(initial_state, mat3_inf_no_imm)[0]
    #inf_no_imm.append(p_inf_no_imm)
    #no_imm.append(p_inf_no_imm)
    # infinite cycles with immigration
    #mat3_inf_with_imm = time_machine.eq_cycle_prob(mat3_with_imm)
    #p_inf_with_imm = np.dot(initial_state, mat3_inf_with_imm)[0]
    #inf_with_imm.append(p_inf_with_imm)
    #with_imm.append(p_inf_with_imm)
    # sequence
    #seqs.append("-".join(seq3))
    labels.append("s-c\nl-p")
    

    print("3")
    print("i-c, s-p")
    (prob5, seq5) = rankings.combine_rankings("rankings/results/eq_seq_ranks_2.txt", "rankings/results/eq_seq_ranks_3.txt", \
                                                "rankings/results/eq_seq_ranks_4.txt", "rankings/results/eq_seq_ranks_5.txt", \
                                                "rankings/results/eq_seq_ranks.txt")
    print(prob5, seq5)
    seq5 = seq5.split(",")
    seq5 = tuple([ab.strip("(").strip().strip(")").strip().strip("'") for ab in seq5])
    # single cycle without immigration
    mat5_no_imm = time_machine.seq_mat(seq5)
    p_single_no_imm = np.dot(initial_state, mat5_no_imm)[0]
    print(p_single_no_imm)
    #single_no_imm.append(p_single_no_imm)
    # single cycle with immigration
    mat5_with_imm = time_machine.seq_mat(seq5, immigration_prob = 0.01)
    p_single_with_imm = np.dot(initial_state, mat5_with_imm)[0]
    #single_with_imm.append(p_single_with_imm)
    # infinite cycles without immigration
    mat5_inf_no_imm = time_machine.eq_cycle_prob(mat5_no_imm)
    p_inf_no_imm = np.dot(initial_state, mat5_inf_no_imm)[0]
    #inf_no_imm.append(p_inf_no_imm)
    #no_imm.append(p_inf_no_imm)
    no_imm.append(prob5)

    # infinite cycles with immigration
    mat5_inf_with_imm = time_machine.eq_cycle_prob(mat5_with_imm)
    p_inf_with_imm = np.dot(initial_state, mat5_inf_with_imm)[0]
    #inf_with_imm.append(p_inf_with_imm)
    with_imm.append(p_inf_with_imm)
    # sequence
    #seqs.append("-".join(seq5))
    labels.append("i-c\ns-p")

    print("4")
    print("i-c, l-p")
    (prob4, seq4) = rankings.combine_rankings("rankings/results/eq_long_period_ranks_2.txt", "rankings/results/eq_long_period_ranks_3.txt",\
                                                "rankings/results/eq_long_period_ranks_4.txt", "rankings/results/eq_long_period_ranks_5.txt", \
                                                "rankings/results/eq_long_period_ranks.txt")
    print(prob4, seq4)
    seq4 = seq4.split(",")
    seq4 = tuple([ab.strip("(").strip().strip(")").strip().strip("'") for ab in seq4])
    # single cycle without immigration
    mat4_no_imm = time_machine.seq_mat(seq4)
    p_single_no_imm = np.dot(initial_state, mat4_no_imm)[0]
    print(p_single_no_imm)
    #single_no_imm.append(p_single_no_imm)
    # single cycle with immigration
    mat4_with_imm = time_machine.seq_mat(seq4, immigration_prob = 0.01)
    p_single_with_imm = np.dot(initial_state, mat4_with_imm)[0]
    #single_with_imm.append(p_single_with_imm)
    # infinite cycles without immigration
    mat4_inf_no_imm = time_machine.eq_cycle_prob(mat4_no_imm)
    p_inf_no_imm = np.dot(initial_state, mat4_inf_no_imm)[0]
    #inf_no_imm.append(p_inf_no_imm)
    #no_imm.append(p_inf_no_imm)
    no_imm.append(prob4)
    # infinite cycles with immigration
    mat4_inf_with_imm = time_machine.eq_cycle_prob(mat4_with_imm)
    p_inf_with_imm = np.dot(initial_state, mat4_inf_with_imm)[0]
    #inf_with_imm.append(p_inf_with_imm)
    with_imm.append(p_inf_with_imm)
    # sequence
    #seqs.append("-".join(seq1))
    labels.append("i-c\nl-p")


    print("starting to plot")
    N = 4
    ind = np.arange(N)
    width = 0.45
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind+.05, no_imm, width, color = "black")
    rects2 = ax.bar(ind+.05+width, with_imm, width, color = "gray")
    #rects3 = ax.bar(ind+.06+2*width, inf_no_imm, width, color=next(color))
    #rects4 = ax.bar(ind+.06+3*width, inf_with_imm, width, color=next(color))

    ax.set_ylabel('final fraction susceptible', fontname = "serif", fontsize = 9)
    ax.set_xticks(ind+.05 + width)
    ax.tick_params(labelsize = 9)
    x_label_list = []
    #count = 1
    #for seq in seqs:
    #    x_label_list.append("{}\n({})".format(seq, count))
    #    count += 1
    ax.set_xticklabels(labels, fontname = "serif", fontsize = 9)

    key = ("no immigration", "immigration")
    ax.legend((rects1[0], rects2[0]), key, frameon = False, prop={"family":"serif", "size": 9})

    autolabel(rects1)
    autolabel(rects2)

    fig = plt.gcf()
    fig.set_size_inches(3.54, 3.54, forward = True)
    plt.savefig("figures/figure_2/strategy_summary.pdf")
    
    '''

    ## Figure 3 (and Table 5) ##
    '''
    file_list = []
    cycle_dict = {}
    #g = open("figures/figure_3/random_initial_states.txt", "r")
    g = open("figures/figure_3/random_initial_states_1000.txt", "r")
    #g = open("figure3/random_initial_states.txt", "r")
    #g = open("random_initial_states_1000.txt", "r")
    i = 0

    
    for line in g:
        i += 1
        print(i)
        
        
        l = line.strip().split(",")
        initial_state = np.array(l).astype("float")
        #print(initial_state)
        
        adaptive_seq, final_geno_dist = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 300, immigration_prob = 0, number = i)
            
        pattern = "-".join(pattern_finder(adaptive_seq, min_reps = 15))
        print("pattern", pattern)
        

        #if pattern == ["FEP", "SAM"] or pattern == ["SAM", "FEP"]:
        #    print(i)
        #    g.close()
        #    break
            
        #break
            
        match_found = False
            
        for cycle in file_list:
            if is_pattern_eq(pattern.split("-"), cycle.split("-")):    
                f = open("figures/figure_3/cycles/{}.txt".format(cycle), "a")
                l = []
                for n in range(len(final_geno_dist)):
                    l.append(str(final_geno_dist[n]))
                f.write(",".join(l) + "\n")
                f.close()
                match_found = True
                break
        if not match_found:
            f = open("figures/figure_3/cycles/{}.txt".format(pattern), "w")
            l = []
            for n in range(len(final_geno_dist)):
                l.append(str(final_geno_dist[n]))
            f.write(",".join(l) + "\n")
            f.close()
            file_list.append(pattern)
    
    
        f = open("figures/figure_4/cycles/AMC-CPR.txt", "r")    
        d = {}
        k = 1
        for line in f:
            #print(k)
            l = line.strip().split(",")
            for i in range(len(l)):
                if i not in d:
                    d[i] = []
                #print(float(l[i]))
                d[i].append(float(l[i]))
            k += 1
        f.close()
        

        mean_list = []
        for j in range(16):
            #space = input("...")
            #print(d[i])
            mean = str(sum(d[j]) / len(d[j]))
            #print(mean, space)
            mean_list.append(mean)
        
        print(mean_list)

        g = open("figures/figure_4/cycles/AMC-CPR-mean.txt", "w")
        g.write(",".join(mean_list))
        g.close()

    
        #len_pattern = len(pattern)
        #pattern = "AMC-CPR"
        
        save = "figures/figure_3/cycles/adaptive_seq_{}.pdf".format(i) 

        # this returns the average growth rate over the last 100 steps
        avg_growth_rate = time_machine.plot_time_series(adaptive_seq, time_machine.cpm, 1, initial_state, savefig = None, immigration_prob = 0)
        #print("avg growth rate:", avg_growth_rate)

        in_dict = False
        for cycle in cycle_dict:
            if is_pattern_eq(pattern.split("-"), cycle.split("-")):
                cycle_dict[cycle]["count"] += 1
                cycle_dict[cycle]["growth rates"].append(avg_growth_rate * 1000)
                # add the new distribution as a new row--each column corresponds to a genotype
                cycle_dict[cycle]["initial geno dist"] = np.append(cycle_dict[cycle]["initial geno dist"], np.array([initial_state]), axis = 0)
                #cycle_dict[cycle]["initial geno dist"] = np.vstack(cycle_dict[cycle]["initial geno dist"], np.array([initial_state]))
                cycle_dict[cycle]["final geno dist"] = np.append(cycle_dict[cycle]["final geno dist"], np.array([final_geno_dist]), axis = 0)
                #cycle_dict[cycle]["final geno"]
                in_dict = True
                #print(cycle_dict[cycle]["final geno dist"].shape)
                #print(cycle_dict[cycle]["final geno dist"]
                print("growth rates:", cycle_dict[cycle]["growth rates"])     
        if not in_dict:
            cycle_dict[pattern] = {"count": 1, "growth rates": [avg_growth_rate * 1000], "random ic growth rates": [], \
                                   "initial geno dist": np.array([initial_state]), "final geno dist": np.array([final_geno_dist])}
            print("growth rates:", cycle_dict[pattern]["growth rates"]) 
               
    g.close()

    num_cycles = len(cycle_dict)
    h = []
    heapq.heapify(h)
    for cycle in cycle_dict:
        print("cycle:", cycle)
        count = cycle_dict[cycle]["count"]
        heapq.heappush(h, (-count, cycle))
        
        # test robustness of cycle 
        for n in range(100):
            print("n:", n)
            if cycle == "None":
                break
            i_s = []
            for m in range(16):
                i_s.append(random.random())
            i_s = np.array(i_s)
            i_s = i_s / sum(i_s)
            # the average growth rate over the last hundred steps
            avg_gr = time_machine.plot_time_series(cycle.split("-"), time_machine.cpm, 300, i_s, savefig = None, immigration_prob = 0)
            cycle_dict[cycle]["random ic growth rates"].append(avg_gr * 1000)

    print("starting to write file 1")
    f1 = open("figures/figure_3/adaptive_cycling_table.csv", "w")
    #f1 = open("figure3/adaptive_cycling_table.csv", "w")
    f1.write("cycle,counts,mean growth rate,max growth rate,mean growth rate (random initial condition),95th percentile growth rate (random initial condition)\n")
    sorted_cycle_list = []
    while h != []:
        print("hi")
        # adaptive cycling table
        neg_count, cycle = heapq.heappop(h)
        rates = cycle_dict[cycle]["growth rates"]
        print("rates:", rates)
        mean_gr = round(sum(rates) / len(rates), 2)
        print("mean_gr:", mean_gr)
        max_gr = round(max(rates), 2)
        print("max_gr:", max_gr)
        if cycle == "None":
            mean_gr_rand_ic = "NA"
            percentile_rand_ic = "NA"
        else:
            rates_rand_ic = cycle_dict[cycle]["random ic growth rates"]
            rates_rand_ic.sort()
            print("rates_rand_ic:", rates_rand_ic)
            mean_gr_rand_ic = round(sum(rates_rand_ic) / len(rates_rand_ic), 2)
            print("mean_gr_rand_ic:", mean_gr_rand_ic)
            percentile_rand_ic = round(rates_rand_ic[94], 2)
            print("percentile_rand_ic:", mean_gr_rand_ic) 
        f1.write("{},{},{},{},{},{}\n".format(cycle, -neg_count, mean_gr, max_gr, mean_gr_rand_ic, percentile_rand_ic))
        sorted_cycle_list.append(cycle)
    f1.close()

    print("starting to write file 2")
    # initial genotype heatmap table
    f2 = open("figures/figure_3/initial_genotype_heatmap_table.txt", "w")
    #f2 = open("figure3/initial_genotype_heatmap_table.txt", "w")
    str_geno = [str(geno) for geno in time_machine.GENO]
    f2.write("," + ",".join(time_machine.GENO) + "\n")
    for cycle in cycle_dict:
        means = np.mean(cycle_dict[cycle]["initial geno dist"], axis = 0).astype("str")
        f2.write(cycle + "," + ",".join(means) + "\n")
    f2.close()


    print("starting to write file 3")
    # final genotype heatmap table 
    f3 = open("figures/figure_3/final_genotype_heatmap_table.txt", "w")
    #f3 = open("figure3/final_genotype_heatmap_table.txt", "w")
    str_geno = [str(geno) for geno in time_machine.GENO]
    f3.write("," + ",".join(str_geno) + "\n")
    for cycle in cycle_dict:
        #for a given cycle, take the means of each of the columns - the mean growth rate for all of the genotypes
        #means = np.mean(cycle_dict[cycle]["final geno dist"], axis = 0).astype("str")
        means = np.mean(cycle_dict[cycle]['final geno dist'], axis = 0).astype("str")
        print(means)
        f3.write(cycle + "," + ",".join(means) + "\n")
    f3.close()
    
    
    for i in ["initial", "final"]:
        #f = open("figure3/{}_genotype_heatmap_table.txt".format(i), "r")
        f = open("figures/figure_3/{}_genotype_heatmap_table.txt".format(i), "r")
        col_labels = f.readline().strip().split(',')[1:]
        num_cycles = 0
        data_list = []
        for line in f:
            num_cycles += 1
            l = line.strip().split(',')
            data_list += [float(num) for num in l[1:]]
        row_labels = list(range(1, num_cycles + 1))
        data = np.array(data_list).reshape((len(row_labels), len(col_labels)))
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                num = round(data[y,x], 3)
                plt.text(x + 0.5, y + 0.5, '{}'.format(num),
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize = 4, fontname = "serif")

        ax.set_ylim(0, num_cycles)
        ax.set_xticks(np.arange(data.shape[1]) + .5, minor = False)
        ax.set_yticks(np.arange(data.shape[0]) + .5, minor = False)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.set_xticklabels(col_labels, minor = False, fontsize = 9, fontname = "serif", rotation = 'vertical')
        ax.set_yticklabels(row_labels, minor = False, fontsize = 9, fontname = "serif")
        
        #cax = ax.imshow(data, interpolation='nearest', cmap=plt.cm.Blues)
        #cbar = fig.colorbar(cax, ticks=[0, 0.5, 1], orientation = "horizontal")
        #bar.ax.set_xticklabels(['0.0', '0.5', '1.0']) 
        
        #fig.tight_layout()
        fig.set_size_inches(3.54, 5, forward = True)
        #plt.savefig("figures/figure_3/{}_genotype_heatmap_text.pdf".format(i))
        plt.savefig("figures/figure_3/{}_genotype_heatmap_text.pdf".format(i))
        plt.close()
        f.close()
        
    '''
    
    ## Figure 4 (and Table 6) ##

    '''
    file_list = []
    cycle_dict = {}
    #g = open("figures/figure_3/random_initial_states.txt", "r")
    g = open("figures/figure_3/random_initial_states_1000.txt", "r")
    #g = open("figure3/random_initial_states.txt", "r")
    #g = open("random_initial_states_1000.txt", "r")
    i = 0

    
    for line in g:
        i += 1
        print(i)
        
        l = line.strip().split(",")
        initial_state = np.array(l).astype("float")
        #print(initial_state)
        
        adaptive_seq, final_geno_dist = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 300, immigration_prob = 0.01, number = i)
            
        pattern = "-".join(pattern_finder(adaptive_seq, min_reps = 15))
        print("pattern", pattern)
        #if pattern == ["FEP", "SAM"] or pattern == ["SAM", "FEP"]:
        #    print(i)
        #    g.close()
        #    break
            
        #break
         
        match_found = False
            
        for cycle in file_list:
            if is_pattern_eq(pattern.split("-"), cycle.split("-")):    
                f = open("figures/figure_4/cycles/{}.txt".format(cycle), "a")
                l = []
                for n in range(len(final_geno_dist)):
                    l.append(str(final_geno_dist[n]))
                f.write(",".join(l) + "\n")
                f.close()
                match_found = True
                break
        if not match_found:
            f = open("figures/figure_4/cycles/{}.txt".format(pattern), "w")
            l = []
            for n in range(len(final_geno_dist)):
                l.append(str(final_geno_dist[n]))
            f.write(",".join(l) + "\n")
            f.close()
            file_list.append(pattern)

        '''
        
    cycle_dict = {}
    #g = open("figures/figure_3/random_initial_states_imm.txt", "r")
    
    g = open("random_initial_states_1000.txt", "r")
    #g = open("figures/figure_3/random_initial_states_1000.txt", "r")

    i = 0
    for line in g:
        i += 1
        print(i)
        l = line.strip().split(",")
        initial_state = np.array(l).astype("float")
    
        adaptive_seq, final_geno_dist = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 300, immigration_prob = 0.01)
        pattern = pattern_finder(adaptive_seq, min_reps = 15)
        len_pattern = len(pattern)
        pattern = "-".join(pattern)
        if pattern == "None":
            save = "figures/figure_4/adaptive_seq_none_{}.pdf".format(i)
        else:
            save = None
        #geno_dist = time_machine.plot_cycle_probs(adaptive_seq, initial_state, time_machine.cpm, savefig = save, immigration_prob = 0.01)
        #final_geno_dist = np.mean(geno_dist[:,-len_pattern:], axis = 1)
        # the average growth rate over the last 100 steps:
        avg_growth_rate = time_machine.plot_time_series(adaptive_seq, time_machine.cpm, 1, initial_state, savefig = None, immigration_prob = 0.01)
        
        in_dict = False
        for cycle in cycle_dict:
            if is_pattern_eq(pattern.split("-"), cycle.split("-")):
                cycle_dict[cycle]["count"] += 1
                cycle_dict[cycle]["growth rates"].append(avg_growth_rate * 1000)
                # add the new distribution as a new row--each column corresponds to a genotype
                cycle_dict[cycle]["initial geno dist"] = np.append(cycle_dict[cycle]["initial geno dist"], np.array([initial_state]), axis = 0)
                cycle_dict[cycle]["final geno dist"] = np.append(cycle_dict[cycle]["final geno dist"], np.array([final_geno_dist]), axis = 0)
                in_dict = True
        if not in_dict:
            cycle_dict[pattern] = {"count": 1, "growth rates": [avg_growth_rate * 1000], "random ic growth rates": [], \
                                   "initial geno dist": np.array([initial_state]), "final geno dist": np.array([final_geno_dist])}
    
    g.close()
    
    h = []
    heapq.heapify(h)
    for cycle in cycle_dict:
        print("cycle:", cycle)
        count = cycle_dict[cycle]["count"]
        heapq.heappush(h, (-count, cycle))
        
        # test robustness of cycle 
        for n in range(100):
            print("n:", n)
            if cycle == "None":
                break
            i_s = []
            for m in range(16):
                i_s.append(random.random())
            i_s = np.array(i_s)
            i_s = i_s / sum(i_s)
            avg_gr = time_machine.plot_time_series(cycle.split("-"), time_machine.cpm, 300, i_s, savefig = None, immigration_prob = 0.01)
            cycle_dict[cycle]["random ic growth rates"].append(avg_gr * 1000)

    print("starting to write file 1")
    f1 = open("figures/figure_4/adaptive_cycling_imm_table.csv", "w")
    #f1 = open("figures/figures_4/adaptive_cycling_imm_table.csv", "w")
    f1.write("cycle,counts,mean growth rate,max growth rate,mean growth rate (random initial condition),95th percentile growth rate (random initial condition)\n")
    sorted_cycle_list = []
    while h != []:
        # adaptive cycling table
        neg_count, cycle = heapq.heappop(h)
        rates = cycle_dict[cycle]["growth rates"]
        mean_gr = round(sum(rates) / len(rates), 2)
        max_gr = round(max(rates), 2)
        if cycle == "None":
            mean_gr_rand_ic = "NA"
            percentile_rand_ic = "NA"
        else:
            rates_rand_ic = cycle_dict[cycle]["random ic growth rates"]
            mean_gr_rand_ic = round(sum(rates_rand_ic) / len(rates_rand_ic), 2)
            rates_rand_ic.sort()
            percentile_rand_ic = round(rates_rand_ic[94], 2) 
        f1.write("{},{},{},{},{},{}\n".format(cycle, -neg_count, mean_gr, max_gr, mean_gr_rand_ic, percentile_rand_ic))
        sorted_cycle_list.append(cycle)
    f1.close()

    print("starting to write file 2")
    # initial genotype heatmap table
    f2 = open("figures/figure_4/initial_genotype_heatmap_imm_table.txt", "w")
    str_geno = [str(geno) for geno in time_machine.GENO]
    f2.write("," + ",".join(time_machine.GENO) + "\n")
    for cycle in sorted_cycle_list:
        means = np.mean(cycle_dict[cycle]["initial geno dist"], axis = 0).astype("str")
        f2.write(cycle + "," + ",".join(means) + "\n")
    f2.close()
    

    print("starting to write file 3")
    # final genotype heatmap table 
    f3 = open("figures/figure_4/final_genotype_heatmap_imm_table.txt", "w")
    str_geno = [str(geno) for geno in time_machine.GENO]
    f3.write("," + ",".join(str_geno) + "\n")
    for cycle in sorted_cycle_list:
        #for a given cycle, take the means of each of the columns - the mean growth rate for all of the genotypes
        means = np.mean(cycle_dict[cycle]["final geno dist"], axis = 0).astype("str")
        f3.write(cycle + "," + ",".join(means) + "\n")
    f3.close()
    
    for i in ["initial", "final"]:
        f = open("figures/figure_4/{}_genotype_heatmap_imm_table.txt".format(i), "r")
        col_labels = f.readline().strip().split(',')[1:]
        num_cycles = 0
        data_list = []
        for line in f:
            num_cycles += 1
            l = line.strip().split(',')
            data_list += [float(num) for num in l[1:]]
        print(num_cycles)
        row_labels = list(range(1, num_cycles + 1))
        data = np.array(data_list).reshape((len(row_labels), len(col_labels)))
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                num = round(data[y,x], 3)
                plt.text(x + 0.5, y + 0.5, '{}'.format(num),
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize = 4, fontname = "serif")

        #print(data.shape[0], data.shape[1])
        ax.set_ylim(0, num_cycles)
        ax.set_xticks(np.arange(data.shape[1]) + .5, minor = False)
        ax.set_yticks(np.arange(data.shape[0]) + .5, minor = False)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.set_xticklabels(scol_labels, minor = False, fontsize = 9, fontname = "serif", rotation = 'vertical')
        ax.set_yticklabels(row_labels, minor = False, fontsize = 9, fontname = "serif")
        
        #cax = ax.imshow(data, interpolation='nearest', cmap=plt.cm.Blues)
        #cbar = fig.colorbar(cax, ticks=[0, 0.5, 1], orientation = "horizontal")
        #bar.ax.set_xticklabels(['0.0', '0.5', '1.0']) 
        
        fig.tight_layout()
        fig.set_size_inches(3.54, 4.5, forward = True)
        plt.savefig("figures/figure_4/{}_genotype_heatmap_imm_text.pdf".format(i))
        plt.close()
        f.close()

    
    
    ## Figure 5 ##


    '''
    f = open("figures/figure_5/table_5.csv", "r")
    f.readline()
    rank_heap = []
    heapq.heapify(rank_heap)
    cycle_dict = {}
    rank_1 = 1
    for line in f:
        l = line.split(",")
        cycle = l[0].strip()
        cycle_dict[cycle] = {"rank 1": rank_1} 
        random_ic_95_percentile = l[-1].strip()
        heapq.heappush(rank_heap, (random_ic_95_percentile, cycle))
        rank_1 += 1

    rank_2 = 1
    while rank_heap != []:
        rand_perc, cycle = heapq.heappop(rank_heap)
        cycle_dict[cycle]["rank 2"] = rank_2
        rank_2 += 1
    #print(cycle_dict)
    cycle_list = [(cycle_dict[cycle]["rank 1"], cycle_dict[cycle]["rank 2"]) for cycle in cycle_dict]
    
    rank_1_list = [cycle[0] for cycle in cycle_list]
    rank_2_list = [cycle[1] for cycle in cycle_list]

    spearman = scipy.stats.spearmanr(rank_1_list, rank_2_list)
    corr = spearman[0]
    p_val = spearman[1]
    print("corr:", corr)
    print("p-val:", p_val)

    fig, ax = plt.subplots()


    ax.plot(rank_1_list, rank_2_list, 'ro', color = "gray")
    ax.set_xlabel("rank with adaptive selection", fontsize = 9, fontname = "serif")
    ax.set_ylabel("rank under random initial conditions", fontsize = 9, fontname = "serif")
    ax.tick_params(axis="both", which="major", labelsize=9)
    ax.tick_params(axis="both", which="minor", labelsize=9)
    x_max = max(rank_1_list) + 1
    y_max = max(rank_2_list) + 1
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)


    #plt.savefig('plots/ranks_top_50_{}_vs_{}_{}.png'.format(x_label, y_label, k))
    fig.set_size_inches(3.54, 3.54, forward = True)
    plt.savefig("figures/figure_5/adaptive_cycling_rank.pdf")
    plt.close()
    f.close()
    '''
    '''
    #figure 3
    initial_state = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    cycle_dict = {}
    g = open("figures/figure_3/random_initial_states.txt", "r")
    #for i in range(100):
    i = 0
    fig, axarr = plt.subplots(1, 3)
    xlabels = [0, 5, 10, 15, 20, 25, 30]
    plt.setp(axarr, xticks=[0, 1/6*30, 2/6*30, 3/6*30, 4/6*30, 5/6*30, 30], xticklabels=xlabels)
    #plt.setp(axarr, xticklabels = xlabels)
    for line in g:
        i += 1
        print(i)
        l = line.strip().split(",")
        initial_state = np.array(l).astype("float")
    
        adaptive_seq = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 150, immigration_prob = 0)
        pattern = pattern_finder(adaptive_seq, min_reps = 10)
        len_pattern = len(pattern)
        pattern = "-".join(pattern)
        if i == 33 or i == 37 or i == 67:
            if i == 33:
                n = 0
                legend = False
            elif i == 37:
                n = 1
                legend = False
            elif i == 67:
                n = 2
                legend = True
            num_cycles = 30
            plot_info = np.zeros((len(time_machine.GENO), num_cycles))
            M = np.identity(len(time_machine.GENO))
            for j in range(num_cycles):
                transition_mat_list = []
                for ab in adaptive_seq:
                    transition_mat_list.append(time_machine.cpm(ab, immigration_prob = 0))
                M = time_machine.cycle_prob(transition_mat_list)
                M = np.dot(M, time_machine.cpm(adaptive_seq[j], immigration_prob = 0))
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

            #xlabels = adaptive_seq
            #print(xlabels)
            axarr[n].tick_params(axis='both', which='major', labelsize=9)
            #xlabels = list(range(1, num_cycles + 1))
            #axarr[n].set_xticklabels(xlabels, ind + 0.4, fontsize = 9, fontname = "serif")
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
    
    plt.savefig("figures/figure_3/FEP_plots.pdf")
    plt.close()
    '''
    
    '''
    # final analysis
    for ip in [0.01]:
        print(ip)
        if ip == 0:
            fig_name = "figures/figure_3/adaptive_cycling_table.csv"
            new_fig_name = "figures/figure_3/adaptive_cycling_table_new.csv"
        else:
            fig_name = "figures/figure_4/adaptive_cycling_imm_table.csv"
            new_fig_name = "figures/figure_4/adaptive_cycling_imm_table_new.csv"

        perc_dict = {}

        g = open("figures/figure_3/random_initial_states_1000.txt", "r")
        #g = open("figure3/random_initial_states.txt", "r")

        f1 = open(fig_name, "r")
        
        f1.readline()
        mean_gr_dict = {}

        for line in f1:
            l = line.strip().split(',')
            pattern = l[0]
            mean_gr = float(l[2])
            mean_gr_dict[pattern] = mean_gr

        f1.close()

        i = 0
        for line in g:
            i += 1
            print(i)
            l = line.strip().split(",")
            initial_state = np.array(l).astype("float")
            
            # find the "real" adaptive sequence pattern:
            adaptive_seq, final_geno_dist = time_machine.adaptive_cycling(initial_state, time_machine.cpm, 300, immigration_prob = ip)
            pattern = pattern_finder(adaptive_seq, min_reps = 15)

            gr_dist = []
            for j in range(1000):
                rand_adaptive_seq, final_geno_dist = time_machine.random_cycling(initial_state, time_machine.cpm, 300, immigration_prob = ip)
                # this returns the average growth rate over the last 100 steps
                avg_growth_rate = time_machine.plot_time_series(rand_adaptive_seq, time_machine.cpm, 1, initial_state, savefig = None, immigration_prob = ip)
                gr_dist.append(avg_growth_rate * 1000)
            gr_dist.sort()
            for k in range(len(gr_dist)):
                print("mean gr:", mean_gr)
                if mean_gr >= gr_dist[pattern][j]:
                    rank += 1
                else:
                    break
            percentile = (rank / len(gr_dist[pattern])) * 100
            if pattern not in perc_dict:
                perc_dict[pattern] = []
            perc_dict[pattern].append(percentile)
            #rand_growth_rate_dist.append(avg_growth_rate * 1000)
        #rand_growth_rate_dist.sort()
        #print(rand_growth_rate_dist)

        # alter adaptive cyclign table:
        f1 = open(fig_name, "r")
        f2 = open(new_fig_name, "w")

        f1.readline()
        f2.write("cycle,counts,mean growth rate,max growth rate,mean growth rate (random initial condition),95th percentile growth rate (random initial condition)\n")

        for line in f1:
            
            line = line.strip().split(',')
            pattern = line[0]
            mean_gr = float(line[2])
            #rank = 0
            #for k in range(len(gr_dist[pattern])):
            #    print("mean gr:", mean_gr)
            #    if mean_gr >= gr_dist[pattern][j]:
            #        rank += 1
            #    else:
            #        break
            #percentile = (rank / len(gr_dist[pattern])) * 100

            #the percentile is averaged over all of the initial conditions coresponding to the pattern
            percentile = sum(perc_dict[pattern]) / len(perc_dict[pattern])
            new_line = "{},{},{},{},{},{}\n".format(line[0],line[1],"{}({})".format(line[2],round(percentile,2)),line[3],line[4],line[5])
            f2.write(new_line + "\n")
        f1.close()
        f2.close() 
    '''
    



    

    