from collections import defaultdict
import pandas as pd
from pareto_mst_plots import get_df, CATEGORIES, add_count_col
import neuron_density
import pylab
from numpy.ma import masked_invalid
from scipy.stats import pearsonr, spearmanr

def size_correlation(df):
    print "---------------size-alpha correlation----------------"
    df2 = df.drop_duplicates(subset='name')
    corr1 = pearsonr(df2['points'], df2['alpha'])
    corr2 = spearmanr(df2['points'], df2['alpha'])
    print "pearson correlation: " + str(corr1) 
    print "spearman correlation: " + str(corr2)

def category_correlation(df, category):
    alphas = df['alpha']
    for unique_val in df[category].unique():
        bit_vec = []
        for val in df[category]:
            bit_vec.append(int(val == unique_val))
        coef, pval = pearsonr(bit_vec, alphas)
        print unique_val, coef, pval

def categories_correlations(df):
    df2 = df.drop_duplicates(subset='name')
    for category in ['species', 'region', 'cell_type']:
        print category
        print '--------------------------------'
        category_correlation(df2, category)

def basic_stats(df):
    df2 = df.drop_duplicates(subset='name')
    total_trials = df2['trials'].sum()
    total_successes = df2['successes'].sum()
    print total_successes, total_trials
    print "p-value", float(total_successes) / total_trials
    print "neural to centroid ratio", infmean(df2['neural_dist'] / df2['centroid_dist'])
    print "neural to random ratio", infmean(df2['neural_dist'] / df2['random_dist'])

    print "dominate percentage", float(df2['dominates'].sum()) / df2['comparisons'].sum()
    df3 = df2[df2['neuron_type'] != 'axon']
    print "dendrite dominate percentage", float(df3['dominates'].sum()) / df3['comparisons'].sum()

    df4 = df2[df2['neural_dist'] < df2['centroid_dist']]
    print "beats centroid", float(df4['neural_dist'].count()) / float(df2['neural_dist'].count())

def infmean(arr):
    return pylab.mean(masked_invalid(arr))

def metadata(df):
    print "unique neurons"
    print len(df['name'].unique())

    for category in CATEGORIES:
        print "unique " + category
        print len(df[category].unique())
        df2 = df.drop_duplicates(subset=['name', category])
        df2 = add_count_col(df2, category)
        f = open(category + '.txt', 'w')
        with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
            print >> f,  df2
        f.close()

def neuron_type_alphas(df):
    df2 = df.drop_duplicates(subset=['name', 'neuron_type'])
    for neuron_type, group in df2.groupby('neuron_type'):
        print neuron_type, pylab.mean(group['alpha'])

def main():
    fname = 'pareto_mst.csv'
    df = get_df()
    metadata(df)
    neuron_type_alphas(df)
    basic_stats(df)
    #categories_correlations(df)
    size_correlation(df)

if __name__ == '__main__':
    main()
