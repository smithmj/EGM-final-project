import scipy.stats
import pandas as pd
#import matplotlib as mpl 
#mpl.rc('font',family='Times New Roman')
import matplotlib.pyplot as plt 


def prob_plot(k, file1, file2, x_label, y_label):
    df1 = pd.read_csv(file1, sep = '\t')
    probs1 = list(df1[df1.columns[1]])
    df2 = pd.read_csv(file2, sep = '\t')
    probs2 = list(df2[df2.columns[1]])

    pearson = scipy.stats.pearsonr(probs1, probs2)
    corr = pearson[0]
    p_val = pearson[1]

    plt.plot(probs1, probs2, 'ro')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title('r = {}'.format(round(corr, 2)))
    plt.savefig('plots/{}_vs_{}_{}.png'.format(x_label, y_label, k))
    plt.close()
    return corr, p_vals

#def rank_plot(k,file1, file2, x_label, y_label, complete_plot = True):
def rank_plot(file1, file2, x_label, y_label, savefig = None, complete_plot = True):
    df1 = pd.read_csv(file1, sep = '\t')
    df2 = pd.read_csv(file2, sep = '\t')

    if complete_plot:
        size = df1.shape[0]
    else:
        size = 50

    rank1 = []
    rank2 = []
    for i in range(size):
        seq = df1[' Sequence'].iloc[[i]].values[0]
        for j in range(df2.shape[0]):
            if df2[' Sequence'].iloc[[j]].values[0] == seq:
                rank1.append(df1['Rank'].iloc[[i]].values[0])
                rank2.append(df2['Rank'].iloc[[j]].values[0])

    spearman = scipy.stats.spearmanr(rank1, rank2)
    corr = spearman[0]
    p_val = spearman[1]

    plt.plot(rank1, rank2, 'ro')
    plt.xlabel(x_label, fontname = "serif")
    plt.ylabel(y_label, fontname = "serif")
    plt.title('r = {}'.format(round(corr, 2)), fontname = "serif")
    #plt.savefig('plots/ranks_top_50_{}_vs_{}_{}.png'.format(x_label, y_label, k))
    if savefig != None:
        plt.savefig(savefig)
    plt.close()
    return corr, p_val

if __name__ == '__main__':
    print("Infinite vs single cycle ")
    #rank_plot(2, 'rankings/results/seq_ranks_2.txt', 'rankings/results/eq_seq_ranks_2.txt', 'Single Cycle Rank', 'Equilibrium Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/seq_ranks_3.txt', 'rankings/results/eq_seq_ranks_3.txt', 'Single Cycle Rank', 'Equilibrium Rank', complete_plot = False)
    corr, p_val = rank_plot('rankings/results/seq_ranks.txt', 'rankings/results/eq_seq_ranks.txt', 'Single cycle rank\n(Short period)', 'Infinite cycle rank\n(Short period)', complete_plot = False)
    print("corr: {}, p-val: {}".format(corr, p_val))
    print("Infinite vs long period")
    #rank_plot(2, 'rankings/results/eq_seq_ranks_2.txt', 'rankings/results/long_period_ranks_2.txt', 'Equilibrium Rank', 'Long Period Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/eq_seq_ranks_3.txt', 'rankings/results/long_period_ranks_3.txt', 'Equilibrium Rank', 'Long Period Rank', complete_plot = False)
    corr, p_val = rank_plot('rankings/results/eq_seq_ranks.txt', 'rankings/results/long_period_ranks.txt', 'Infinite cycle rank\n(Short period)', 'Single cycle rank\n(Long period)', complete_plot = False)
    print("corr: {}, p-val: {}".format(corr, p_val))
    print("Single vs long period")
    #rank_plot(2, 'rankings/results/seq_ranks_2.txt', 'rankings/results/long_period_ranks_2.txt',  'Single Cycle Rank', 'Long Period Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/seq_ranks_3.txt', 'rankings/results/long_period_ranks_3.txt',  'Single Cycle Rank', 'Long Period Rank', complete_plot = False)
    corr, p_val = rank_plot('rankings/results/seq_ranks.txt', 'rankings/results/long_period_ranks.txt',  'Short period rank\n(Single cycle)', 'Long period rank\n(Single cycle)', complete_plot = False)
    print("corr: {}, p-val: {}".format(corr, p_val))
    #print("4")
    #rank_plot(2, 'rankings/results/eq_seq_ranks_2.txt', 'rankings/results/eq_ranks_5perc_imm_2.txt', 'Equilibrium Rank', 'Equilibrium with Immigration Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/eq_seq_ranks_3.txt', 'rankings/results/eq_ranks_5perc_imm_3.txt', 'Equilibrium Rank', 'Equilibrium with Immigration Rank', complete_plot = False)
    #rank_plot('rankings/results/eq_seq_ranks.txt', 'rankings/results/eq_seq_ranks_5perc_imm.txt', 'Infinite cycle rank', 'Infinite cycle with immigration rank', complete_plot = False)
    #print("5")
    #rank_plot(2, 'rankings/results/seq_ranks_2.txt', 'rankings/results/seq_ranks_5perc_imm_2.txt', 'Single Cycle Rank', 'Single Cycle with Immigration Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/seq_ranks_3.txt', 'rankings/results/seq_ranks_5perc_imm_3.txt', 'Single Cycle Rank', 'Single Cycle with Immigration Rank', complete_plot = False)
    #rank_plot('rankings/results/seq_ranks.txt', 'rankings/results/seq_ranks_5perc_imm.txt', 'Single cycle rank\n(Short period)', 'Single cycle with immigration rank\n(Short period)', complete_plot = False)
    #print("6")
    #rank_plot(2, 'rankings/results/long_period_ranks_2.txt', 'rankings/results/long_period_ranks_5perc_imm_2.txt', 'Long Period Rank', 'Long Period with Immigration Rank', complete_plot = False)
    #rank_plot(3, 'rankings/results/long_period_ranks_3.txt', 'rankings/results/long_period_ranks_5perc_imm_3.txt', 'Long Period Rank', 'Long Period with Immigration Rank', complete_plot = False)
    #rank_plot('rankings/results/long_period_ranks.txt', 'rankings/results/long_period_ranks_5perc_imm.txt', 'Single cycle rank\n(Long period)', 'Single cycle with immigration rank\n(Long period)', complete_plot = False)
    