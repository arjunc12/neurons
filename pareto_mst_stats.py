from collections import defaultdict
import pandas as pd
from pareto_mst_plots import COLUMNS
import pylab
from numpy.ma import masked_invalid
from scipy.stats import pearsonr

def category_correlation(df, category):
    alphas = df['alpha']
    for unique_val in df[category].unique():
        bit_vec = []
        for val in df[category]:
            bit_vec.append(int(val == unique_val))
        coef, pval = pearsonr(bit_vec, alphas)
        print unique_val, coef, pval

def categories_correlations(df):
    for category in ['species', 'region', 'cell_type']:
        print category
        print '--------------------------------'
        category_correlation(df, category)

def basic_stats(df):
    total_trials = df['trials'].sum()
    total_successes = df['successes'].sum()
    print "p-value", float(total_successes) / total_trials
    print "neural to centroid ratio", infmean(df['neural_dist'] / df['centroid_dist'])
    print "neural to random ratio", infmean(df['neural_dist'] / df['random_dist'])
    print "dominate percentage", float(df['dominates'].sum()) / df['comparisons'].sum()

    df2 = df[df['neural_dist'] < df['centroid_dist']]
    print "beats centroid", float(df2['neural_dist'].count()) / float(df['neural_dist'].count())

def infmean(arr):
    return pylab.mean(masked_invalid(arr))

def main():
    fname = 'pareto_mst.csv'
    df = pd.read_csv(fname, names=COLUMNS)
    #basic_stats(df)
    categories_correlations(df)

if __name__ == '__main__':
    main()
