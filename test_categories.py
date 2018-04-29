import pandas as pd
from scipy.stats import ttest_ind, ks_2samp
from pareto_steiner_stats import get_dfs, OUTPUT_FILE, MODELS_FILE, CATEGORIES_FILE
import argparse
from itertools import combinations
import pylab

def test_cat_vals(df, category, vals):
    df2 = df[df[category].isin(vals)]
    df2 = df2.drop_duplicates(subset=['neuron name', 'neuron type'])
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
    parser.add_argument('-l', '--labs', nargs='+', default=None)
    parser.add_argument('-of', '--output_file', default=OUTPUT_FILE)
    parser.add_argument('-mf', '--models_file', default=MODELS_FILE)
    parser.add_argument('-cf', '--categories_file', default=CATEGORIES_FILE)
    parser.add_argument('--synthetic', action='store_true')

    args = parser.parse_args()
    species = args.species
    regions = args.regions
    cell_types = args.cell_types
    labs = args.labs
    output_file = args.output_file
    models_file = args.models_file
    categories_file = args.categories_file
    synthetic = args.synthetic

    if synthetic:
        output_file = output_file.replace('.csv', '_synthetic.csv')
        models_file = models_file.replace('.csv', '_synthetic.csv')
    
    models_df, categories_df =  get_dfs(output_file=OUTPUT_FILE,\
                                        categories_file=CATEGORIES_FILE,\
                                        models_file=MODELS_FILE)

    for category, vals in zip(['species', 'region', 'cell type', 'lab'],\
                              [species, regions, cell_types, labs]):
        if vals != None:
            test_cat_vals(categories_df, category, vals)

if __name__ == '__main__':
    main()
