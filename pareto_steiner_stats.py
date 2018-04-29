from collections import defaultdict
import pandas as pd
import neuron_density
import pylab
from numpy.ma import masked_invalid
from scipy.stats import *
import os
import seaborn as sns
from itertools import combinations
import numpy as np
from scipy.stats import entropy, binom_test, ttest_1samp
from numpy.linalg import norm
import numpy as np
from stats_utils import *
import argparse

FIGS_DIR = 'steiner_stats'

TEST_NEW_FUNCTION = False

OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output'
OUTPUT_FNAME = 'pareto_steiner.csv'
OUTPUT_FILE = '%s/%s' % (OUTPUT_DIR, OUTPUT_FNAME)
MODELS_FNAME = 'models.csv'
MODELS_FILE = '%s/%s' % (OUTPUT_DIR, MODELS_FNAME)

CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES = ['cell type', 'species', 'region', 'neuron type', 'lab']

METADATA_DIR = '/iblsn/data/Arjun/neurons/metadata'

CATEGORY_MIN_COUNTS = {'species': 100, 'cell type' : 500, 'region' : 500, 'neuron type' : 100, 'lab' : 500}
MIN_COUNT = 50

MIN_POINTS = 100
MIN_ALPHA = 0
MAX_ALPHA = 1

MAX_POINTS = float("inf")

LOG_DIST = True

def count_duplicate_rows(df):
    all_rows = len(df.index)
    df2 = df.drop_duplicates()
    unique_rows = len(df2.index)
    return all_rows - unique_rows

def add_count_col(df, categories):
    return df.groupby(categories).size().reset_index(name='count')

def remove_small_counts(df, categories, min_count=MIN_COUNT):
    df2 = add_count_col(df, categories)
    df2 = pd.merge(df, df2)
    df2 = df2[df2['count'] >= min_count]
    return df2

def get_dfs(output_file=OUTPUT_FILE, categories_file=CATEGORIES_FILE,\
            models_file=MODELS_FILE):
    output_df = pd.read_csv(output_file, skipinitialspace=True)
    output_df = output_df[output_df['points'] >= MIN_POINTS]
    output_df = output_df[(output_df['alpha'] >= MIN_ALPHA) & (output_df['alpha'] <= MAX_ALPHA)]

    models_df = pd.read_csv(models_file, skipinitialspace=True)
    models_cols = models_df.columns.values
    models_df = pd.merge(models_df, output_df, on=['neuron name', 'neuron type'])
    models_df = models_df[models_cols]
    
    neural_df = models_df[models_df['model'] == 'neural']
    neural_df = neural_df[['neuron name', 'neuron type', 'dist']]

    cat_df = pd.read_csv(categories_file, skipinitialspace=True)
    categories_df = pd.merge(output_df, cat_df, on='neuron name')
    categories_df = pd.merge(categories_df, neural_df, on=['neuron name', 'neuron type'])

    return models_df, categories_df

def get_filtered_df(df=None):
    if df is None:
        df = get_df()
    filtered_df = df.copy()
    for category in CATEGORIES:
        filtered_df = remove_small_counts(filtered_df, category,\
                                          min_count=CATEGORY_MIN_COUNTS[category])
        filtered_df.drop('count', inplace=True, axis=1)
    return filtered_df

def alpha_counts(df, category, cat_value, alphas=None):
    alpha_values = df['alpha'][df[category] == cat_value]
    alpha_values = pylab.array(alpha_values)
    alpha_values = np.around(alpha_values, decimals=2)
    alpha_values = list(alpha_values)
    alpha_values = map(lambda x : round(x, 2), alpha_values)
    
    if alphas == None:
        delta = 0.01
        alphas = pylab.arange(0, 1 + delta, delta)
        alphas = list(alphas)
        alphas = map(lambda x : round(x, 2), alphas)
        
    counts_dict = defaultdict(int)
    for alpha_value in alpha_values:
        counts_dict[alpha_value] += 1

    counts = []
    for alpha in alphas:
        count = counts_dict[alpha]
        counts.append(count)
    
    return counts

def all_counts(df, category, alphas=None):
    counts = {}
    for cat_val in df[category].unique():
        counts[cat_val] = alpha_counts(df, category, cat_val, alphas)

    return counts

def JSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (entropy(_P, _M) + entropy(_Q, _M))

def pseudo_kld(counts1, counts2):
    assert len(counts1) == len(counts2)
    pseudocounts1 = []
    pseudocounts2 = []
    for i in xrange(len(counts1)):
        c1 = counts1[i]
        c2 = counts2[i]
        '''
        if c1 == 0 or c2 == 0:
            continue
        else:
            pseudocounts1.append(c1)
            pseudocounts2.append(c2)
        '''
        pseudocounts1.append(max(c1, PSEUDOCOUNT))
        pseudocounts2.append(max(c2, PSEUDOCOUNT))

    kld = entropy(pseudocounts1, pseudocounts2)

    return kld

def normalize_distribution(dist):
    return pylab.array(dist, dtype=np.float64) / sum(dist)

def total_variation_distance(dist1, dist2):
    d1 = normalize_distribution(dist1)
    d2 = normalize_distribution(dist2)

    return np.abs(d1 - d2).max()

def hellinger_distance(dist1, dist2):
    assert len(dist1) == len(dist2)

    d1 = normalize_distribution(dist1)
    d2 = normalize_distribution(dist2)
    
    d1 **= 0.5
    d2 **= 0.5

    dist = d1 - d2
    dist = np.dot(dist, dist)
    dist /= 2
    dist **= 0.5

    return dist

def make_dist_frame(df, category, alphas=None, dist_func=pseudo_kld):
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type', category])
    df2 = remove_small_counts(df2, category,\
                              min_count=CATEGORY_MIN_COUNTS[category])

    counts = all_counts(df2, category, alphas)
    
    cat1 = []
    cat2 = []
    dist_vals = []
    for val1, val2 in combinations(counts.keys(), 2):
        counts1 = counts[val1]
        counts2 = counts[val2]
        #kld = entropy(counts1, counts2)
        dist1 = dist_func(counts1, counts2)
        dist2 = dist_func(counts2, counts1)
        #kld = jsdiv(counts1, counts2)
        cat1 += [val1, val2]
        cat2 += [val2, val1]
        dist_vals += [dist1, dist2]

    dist_frame = pd.DataFrame()
    dist_frame[category + '1'] = cat1
    dist_frame[category + '2'] = cat2
    dist_frame['distance'] = dist_vals

    return dist_frame

#DIST_FUNCS = [pseudo_kld, hellinger_distance]
DIST_FUNCS = [hellinger_distance, JSD, total_variation_distance]
DIST_FUNC_NAMES = {pseudo_kld : 'kld', hellinger_distance : 'hellinger',\
                   JSD : 'jsd', total_variation_distance: 'tvd'}

def dist_heat(df, category, alphas=None, dist_func=pseudo_kld, outdir=FIGS_DIR):
    dist_frame = make_dist_frame(df, category, alphas, dist_func)
    dist_frame = dist_frame.pivot(category + '1', category + '2', 'distance')
    pylab.figure()
    ax = sns.heatmap(dist_frame, vmin=0, vmax=1)
    pylab.xticks(rotation='vertical')
    pylab.yticks(rotation='horizontal')
    pylab.tight_layout()
    pylab.savefig('%s/%s_heat_%s.pdf' % (outdir, DIST_FUNC_NAMES[dist_func],\
                                         category.replace(' ', '_')),
                   format='pdf')
    pylab.close()

def kld_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=pseudo_kld)

def hellinger_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=hellinger_distance)

def jsd_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=JSD)

def dist_heats(df, categories, dist_funcs, alphas=None, outdir=FIGS_DIR):
    for category in categories:
        for dist_func in dist_funcs:
            dist_heat(df, category, alphas=alphas, dist_func=dist_func, outdir=outdir)

def alphas_heat(df, categories, outdir=FIGS_DIR):
    for cat1, cat2 in combinations(categories, 2):
        df2 = df.drop_duplicates(subset=['neuron name', 'neuron type', cat1, cat2])
        df2 = remove_small_counts(df2, [cat1, cat2])
        df2 = df2.groupby([cat1, cat2], as_index=False).agg({'alpha' : pylab.mean})
        data = df2.pivot(cat1, cat2, 'alpha')
        pylab.figure()
        ax = sns.heatmap(data, vmin=0, vmax=1)
        pylab.xticks(rotation='vertical')
        pylab.yticks(rotation='horizontal')
        pylab.tight_layout()
        pylab.savefig('%s/%s_%s_alphas_heat.pdf' % (outdir,\
                                                    cat1.replace(' ', '_'),\
                                                    cat2.replace(' ', '_')),\
                      format='pdf')
        pylab.close()

def cat_to_num(categories):
    unique_categories = set()
    index = 1
    cat_map = {}
    cat_nums = []
    for category in categories:
        if category not in unique_categories:
            cat_map[category] = index
            unique_categories.add(category)
            index += 1
        else:
            assert category in cat_map
        cat_nums.append(cat_map[category])
    return cat_nums

def alpha_distribution(df, categories, plot_func, plot_descriptor, outdir=FIGS_DIR):
    for category in categories:
        subset_cols = ['neuron name', 'neuron type', 'alpha']
        if category != 'neuron type':
            subset_cols.append(category)
        df2 = df.drop_duplicates(subset=subset_cols)
        df2 = remove_small_counts(df2, category,\
                                  min_count=CATEGORY_MIN_COUNTS[category])
        cat_vals = []
        medians = []
        for name, group in df2.groupby(category):
            cat_vals.append(name)
            medians.append(pylab.median(group['alpha']))
        
        cat_vals = pylab.array(cat_vals)
        mean = pylab.array(medians)
        order = pylab.argsort(medians)
        order = cat_vals[order]

        pylab.figure()
        sns.set()
        dist_plot = plot_func(x='alpha', y=category, data=df2, orient='h', order=order)
        dist_plot.tick_params(axis='y', labelsize=20)
        #pylab.tight_layout()
        pylab.xlabel('alpha', fontsize=20)
        pylab.ylabel(category, fontsize=20)
        pylab.savefig('%s/%s_alphas_%s.pdf' % (outdir, category.replace(' ', '_'), plot_descriptor),
                       format='pdf', bbox_inches='tight')
        pylab.close()

def cluster_alphas(df, identifiers, outdir=FIGS_DIR):
    alpha_distribution(df, identifiers, sns.stripplot, 'cluster', outdir=outdir)

def boxplot_alphas(df, identifiers, outdir=FIGS_DIR):
    alpha_distribution(df, identifiers, sns.boxplot, 'box', outdir=outdir)

def violin_alphas(df, identifiers, outdir=FIGS_DIR):
    alpha_distribution(df, identifiers, sns.violinplot, 'violin', outdir=outdir)

def swarm_alphas(df, identifiers, outdir=FIGS_DIR):
    alpha_distribution(df, identifiers, sns.swarmplot, 'swarm', outdir=outdir)
    

def category_dists(df, categories, outdir=FIGS_DIR):
    for category in categories:
        df2 = df.drop_duplicates(subset=['neuron name', 'neuron type', category])
        df2 = remove_small_counts(df2, category,\
                                  min_count=CATEGORY_MIN_COUNTS[category])
        cat_vals = []
        cat_means = []
        for cat_val, group in df2.groupby(category):
            cat_vals.append(cat_val)
            cat_mean = pylab.mean(group['dist'])
            cat_means.append(cat_mean)
        order = pylab.argsort(cat_means)
        cat_vals = pylab.array(cat_vals)
        sorted_vals = cat_vals[order]
        pylab.figure()
        sns.set()
        dist_plot = sns.barplot(x=category, y='dist', data=df2, order=sorted_vals)
        pylab.xticks(rotation='vertical')
        pylab.tight_layout()
        pylab.savefig('%s/pareto_dists_%s.pdf' % (outdir,\
                                                  category.replace(' ', '_')),\
                                                  format='pdf')

def scatter_dists(models_df, outdir=FIGS_DIR):
    df = models_df[['neuron name', 'neuron type', 'model', 'dist']]
    model_dists = defaultdict(list)
    for name, group in df.groupby(['neuron name', 'neuron type']):
        for model, group2 in group.groupby('model'):
            model_dists[model].append(pylab.mean(group2['dist']))

    order = pylab.argsort(model_dists['neural'])
    x = pylab.arange(len(model_dists['neural']))
    pylab.figure()
    sns.set()
    pylab.tight_layout()
    for model, dists in model_dists.iteritems():
        dists = pylab.array(dists)
        y = dists[order]
        if LOG_DIST:
            y = pylab.log10(y)
        pylab.scatter(x, y, label=model)
    pylab.xlabel('neuron index', fontsize=25)
    ylab = ''
    if LOG_DIST:
        ylab += 'log-'
    ylab += 'distance to Pareto front'
    pylab.ylabel(ylab, fontsize=25)
    pylab.legend(ncol=2)
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    pylab.savefig('%s/pareto_dists.pdf' % outdir, format='pdf')

def alphas_hist(df, outdir=FIGS_DIR, categories=None):
    subset_cols = ['neuron name', 'neuron type']
    if categories != None:
        for category in categories:
            if category != 'neuron type':
                subset_cols.append(category)
    df2 = df.drop_duplicates(subset_cols)
    
    alphas = None
    weights = None
    labels = None
    
    if categories == None:
        alphas = list(df2['alpha'])
        print "all neurons mean alpha", pylab.mean(alphas)
        weights = pylab.ones_like(alphas) / len(alphas)
    else:
        alphas = []
        weights = []
        labels = []
        for name, group in df2.groupby(categories):
            cat_alphas = group['alpha']
            print name + " neurons mean alpha", pylab.mean(cat_alphas)
            cat_weights = pylab.ones_like(cat_alphas) / len(cat_alphas)
            alphas.append(cat_alphas)
            weights.append(cat_weights)
            labels.append(name)
    
    pylab.figure()
    sns.set()
    if labels == None:
        pylab.hist(alphas, range=(0, 1), weights=weights)
    else:
        pylab.hist(alphas, range=(0, 1), weights=weights, label=labels)
        pylab.legend()
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    pylab.xlabel('alpha')
    pylab.ylabel('proportion')
    pylab.tight_layout()
    name = 'alphas_hist'
    if categories != None:
        cat_str = '_'.join(categories)
        cat_str = cat_str.replace(' ', '_')
        name += '_' + cat_str

    outname = '%s/%s.pdf' % (outdir, name)
    outname = outname.replace(' ', '_')
    pylab.savefig('%s/%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

def neuron_types_hist(df, outdir=FIGS_DIR):
    alphas_hist(df, outdir, ['neuron type'])

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

def null_models_analysis(models_df):
    print "-----------------------------------------------------"
    df2 = models_df[models_df['model'] != 'neural']
    for model, group in df2.groupby('model'):
        print '***%s***' % model
        ntrials = len(group['success'])
        success_rate = pylab.mean(group['success'])
        print "success rate", success_rate, "trials", len(group['success'])
        print "binomial p-value", binom_test(success_rate, ntrials)
        print "neural to %s ratio" % model, pylab.mean(group['ratio'])
        print "t-test p-value", ttest_1samp(group['ratio'], popmean=1)

def infmean(arr):
    return pylab.mean(masked_invalid(arr))

def metadata(df):
    print "-----------------------------------------------------"
    print "unique neurons"
    print len(df.groupby(['neuron name', 'neuron type']))

    for category in CATEGORIES:
        print "unique " + category
        print len(df[category].unique())
        df2 = df.drop_duplicates(subset=['neuron name', category])
        df2 = add_count_col(df2, category)
        df2 = df2[df2['count'] >= 25]
        category_str = category.replace(' ', '_')
        f = open('%s/%s.txt' % (METADATA_DIR, category_str), 'w')
        with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
            print >> f,  df2
        f.close()

def neuron_type_alphas(df):
    print "-----------------------------------------------------"
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type'])
    types = []
    dists = []
    for neuron_type, group in df2.groupby('neuron type'):
        print neuron_type, pylab.mean(group['alpha']), '+/-', pylab.std(group['alpha'], ddof=1)
        types.append(neuron_type)
        dists.append(pylab.array(group['alpha']))
    indices = range(len(types))
    for idx1, idx2 in combinations(indices, 2):
        type1, type2 = types[idx1], types[idx2]
        dist1, dist2 = dists[idx1], dists[idx2]
        print type1 + ' vs. ' + type2
        #print ttest_ind(dist1, dist2, equal_var=False)
        #print mannwhitneyu(dist1, dist2, alternative='two-sided')
        print ks_2samp(dist1, dist2)

def vals_correlation(df, val1, val2, outdir=FIGS_DIR, xtransform=None,\
                     ytransform=None, logtransform=False):
    print "-----------------------------------------------------"
    print "%s-%s correlation" % (val1, val2)
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type'])
    v1 = df2[val1]
    v2 = df2[val2]
    
    if logtransform:
        xtransform = pylab.log10
        ytransform = pylab.log10
    
    if xtransform != None:
        v1 = xtransform(v1)
    if ytransform != None:
        v2 = ytransform(v2)

    print pearsonr(v1, v2)
    print spearmanr(v1, v2)

    regression_df = df2.copy()
    add_regression_cols(regression_df, val1, val2, xtransform=xtransform,\
                        ytransform=ytransform)

    sns.set()
    pylab.figure()
    nrows = len(df2['neuron type'].unique()) + 1
    pylab.subplot(nrows, 1, 1)
    pylab.scatter(v1, v2)
    x = v1
    y = regression_df['%s_hat' % val2]
    order = pylab.argsort(x)
    x = x[order]
    y = y[order]
    pylab.plot(x, y, c='g')

    row = 2
    for neuron_type, group in df2.groupby('neuron type'):
        print neuron_type
        pylab.subplot(nrows, 1, row)
        row += 1
        v1 = pylab.array(group[val1])
        v2 = pylab.array(group[val2])
        if xtransform != None:
            v1 = xtransform(v1)
        if ytransform != None:
            v2 = ytransform(v2)
        print pearsonr(v1, v2)
        print spearmanr(v1, v2)
    
        regression_df = group.copy()
        add_regression_cols(regression_df, val1, val2, xtransform=xtransform,\
                            ytransform=ytransform)
        
        pylab.scatter(v1, v2)
        x = pylab.array(v1)
        y = pylab.array(regression_df['%s_hat' % val2])
        order = pylab.argsort(x)
        x = x[order]
        y = y[order]
        pylab.plot(x, y, c='g')
    
    pylab.tight_layout()
    figname = '%s/%s_%s.pdf' % (outdir, val1, val2)
    pylab.savefig('%s/%s_%s.pdf' % (outdir, val1, val2), format='pdf')

    pylab.close()

def size_dist_correlation(df, outdir=FIGS_DIR):
    vals_correlation(df, 'points', 'dist', outdir=outdir, logtransform=True)

def alpha_dist_correlation(df, outdir=FIGS_DIR):
    vals_correlation(df, 'alpha', 'dist', outdir=outdir, logtransform=True)

def size_alpha_correlation(df, outdir=FIGS_DIR):
    vals_correlation(df, 'alpha', 'points', outdir=outdir, logtransform=True)

def truncation_hist(df, outdir=FIGS_DIR):
    df2 = df[df['neuron type'].isin(['axon', 'truncated axon'])]
    alphas = []
    weights = []
    labels = []
    for neuron_type, group in df2.groupby('neuron type'):
        alpha = group['alpha']
        weight = pylab.ones_like(alpha) / len(alpha)
        alphas.append(alpha)
        weights.append(weight)
        labels.append(neuron_type)

    pylab.figure()
    sns.set()
    
    pylab.hist(alphas, range=(0, 1), weights=weights, label=labels)
    pylab.legend()
    
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    
    pylab.xlabel('alpha')
    pylab.ylabel('proportion')
    pylab.tight_layout()
    
    name = 'truncation_hist'
    outname = '%s/%s.pdf' % (outdir, name)
    outname = outname.replace(' ', '_')
    pylab.savefig('%s/%s.pdf' % (outdir, name), format='pdf')
    pylab.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-od', '--output_dir', default=OUTPUT_DIR)
    parser.add_argument('-o', '--output_fname', default=OUTPUT_FNAME)
    parser.add_argument('-m', '--models_fname', default=MODELS_FNAME)
    parser.add_argument('-c', '--categories_file', default=CATEGORIES_FILE)
    parser.add_argument('-f', '--figs_dir', default=FIGS_DIR)
    parser.add_argument('--synthetic', action='store_true')
    args = parser.parse_args()
    
    output_dir = args.output_dir
    output_fname = args.output_fname
    models_fname = args.models_fname
    categories_file = args.categories_file 
    figs_dir = args.figs_dir
    synthetic = args.synthetic

    if synthetic:
        output_fname = output_fname.replace('.csv', '_synthetic.csv')
        models_fname = models_fname.replace('.csv', '_synthetic.csv')
        figs_dir += '_synthetic'

    output_file = '%s/%s' % (output_dir, output_fname)
    models_file = '%s/%s' % (output_dir, models_fname)
    
    models_df, categories_df = get_dfs(output_file=output_file,\
                                       categories_file=categories_file,\
                                       models_file=models_file)

    if TEST_NEW_FUNCTION:
        
        return None
    
    metadata(categories_df)
    print "-----------------------------------------------------"
    print "mean neural dist", pylab.mean(models_df['dist'][models_df['model'] == 'neural'])
    null_models_analysis(models_df)
    neuron_type_alphas(categories_df)
    os.system('mkdir -p %s' % figs_dir)
    scatter_dists(models_df, outdir=figs_dir)
    boxplot_alphas(categories_df, CATEGORIES, outdir=figs_dir) 
    alphas_hist(categories_df, outdir=figs_dir)
    neuron_types_hist(categories_df, outdir=figs_dir)
    alphas_heat(categories_df, CATEGORIES, outdir=figs_dir)
    dist_heats(categories_df, CATEGORIES, DIST_FUNCS, outdir=figs_dir)
    category_dists(categories_df, CATEGORIES, outdir=figs_dir)
    truncation_hist(categories_df, outdir=figs_dir)
    size_dist_correlation(categories_df, outdir=figs_dir)
    alpha_dist_correlation(categories_df, outdir=figs_dir)
    size_alpha_correlation(categories_df, outdir=figs_dir)
    
if __name__ == '__main__':
    main()
