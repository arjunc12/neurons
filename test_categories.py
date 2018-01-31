import pandas as pd
from scipy.stats import ttest_ind, ks_2samp
from pareto_steiner_stats import get_df
import argparse
from itertools import combinations
import pylab

def test_cat_vals(df, category, vals):
    df2 = df.drop_duplicates(subset='name')
    df2 = df2[df2[category].isin(vals)]
    alphas = df2['alpha']
    for val1, val2 in combinations(vals, 2):
        sample1 = alphas[df2[category] == val1]
        sample2 = alphas[df2[category] == val2]
        degf = len(sample1) + len(sample2) - 2
        print val1, val2, len(sample1), len(sample2), degf
        print "t test", ttest_ind(sample1, sample2, equal_var=False)
        print "ks test", ks_2samp(sample1, sample2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', nargs='+', default=None)
    parser.add_argument('-r', '--regions', nargs='+', default=None)
    parser.add_argument('-c', '--cell_types', nargs='+', default=None)
    parser.add_argument('-d', '--data_file', default='pareto_steiner.csv')

    args = parser.parse_args()
    species = args.species
    regions = args.regions
    cell_types = args.cell_types
    data_file = args.data_file

    df = get_df(data_file)

    for category, vals in zip(['species', 'region', 'cell_type'], [species, regions, cell_types]):
        if vals != None:
            test_cat_vals(df, category, vals)

if __name__ == '__main__':
    main()
