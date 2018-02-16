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
from scipy.stats import entropy
from numpy.linalg import norm
import numpy as np
from stats_utils import *
import argparse

OUTDIR = 'steiner_stats'

TEST_NEW_FUNCTION = True

OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output'
OUTPUT_FILE = '%s/pareto_steiner.csv' % OUTPUT_DIR
MODELS_FILE = '%s/models.csv' % OUTPUT_DIR

CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES = ['cell type', 'species', 'region', 'neuron type', 'lab']

MIN_COUNT = 50

MIN_POINTS = 100
MAX_POINTS = float("inf")

LOG_DIST = True

def add_count_col(df, categories):
    return df.groupby(categories).size().reset_index(name='count')

def remove_small_counts(df, categories):
    df2 = add_count_col(df, categories)
    df2 = pd.merge(df, df2)
    df2 = df2[df2['count'] >= MIN_COUNT]
    return df2

def get_dfs(output_file=OUTPUT_FILE, categories_file=CATEGORIES_FILE,\
            models_file=MODELS_FILE):
    df1 = pd.read_csv(output_file, skipinitialspace=True)
    df2 = pd.read_csv(categories_file, skipinitialspace=True)
    categories_df = pd.merge(df1, df2)

    models_df = pd.read_csv(models_file, skipinitialspace=True)

    return categories_df, models_df

def get_filtered_df(df=None):
    if df is None:
        df = get_df()
    filtered_df = df.copy()
    for category in CATEGORIES:
        filtered_df = remove_small_counts(filtered_df, category)
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
    df2 = df.drop_duplicates(subset=['name', category])
    df2 = remove_small_counts(df2, category)

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

def dist_heat(df, category, alphas=None, dist_func=pseudo_kld, outdir=OUTDIR):
    dist_frame = make_dist_frame(df, category, alphas, dist_func)
    dist_frame = dist_frame.pivot(category + '1', category + '2', 'distance')
    pylab.figure()
    ax = sns.heatmap(dist_frame, vmin=0, vmax=1)
    ax.tick_params(axis='x', rotation=90)#, labelsize=20)
    ax.tick_params(axis='y', rotation=0)#, labelsize=20)
    pylab.tight_layout()
    pylab.savefig('%s/%s_heat_%s.pdf' % (outdir, DIST_FUNC_NAMES[dist_func], category),
                   format='pdf')#, bbox_inches='tight')
    pylab.close()

def kld_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=pseudo_kld)

def hellinger_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=hellinger_distance)

def jsd_heat(df, category, alphas=None):
    return dist_heat(df, category, alphas=alphas, dist_func=JSD)

def dist_heats(df, categories, dist_funcs, alphas=None, outdir=OUTDIR):
    for category in categories:
        for dist_func in dist_funcs:
            dist_heat(df, category, alphas=alphas, dist_func=dist_func, outdir=outdir)

def alphas_heat(df, categories, outdir=OUTDIR):
    for cat1, cat2 in combinations(categories, 2):
        df2 = df.drop_duplicates(subset=['name', cat1, cat2])
        df2 = remove_small_counts(df2, [cat1, cat2])
        df2 = df2.groupby([cat1, cat2], as_index=False).agg({'alpha' : pylab.mean})
        data = df2.pivot(cat1, cat2, 'alpha')
        pylab.figure()
        ax = sns.heatmap(data, vmin=0, vmax=1)
        pylab.xticks(rotation=90)#, size=20)
        pylab.tight_layout()
        pylab.savefig('%s/%s_%s_alphas_heat.pdf' % (outdir, cat1, cat2),\
                      format='pdf')#, bbox_inches='tight')
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

def alpha_distribution(df, categories, plot_func, plot_descriptor, outdir=OUTDIR):
    for category in categories:
        df2 = df.drop_duplicates(subset=['name', category])
        df2 = remove_small_counts(df2, category)
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
        dist_plot = plot_func(x='alpha', y=category, data=df2, orient='h', order=order)
        #dist_plot.tick_params(axis='y', labelsize=20)
        pylab.tight_layout()
        pylab.savefig('%s/%s_alphas_%s.pdf' % (outdir, category, plot_descriptor),
                       format='pdf')#, bbox_inches='tight')
        pylab.close()

def cluster_alphas(df, identifiers, outdir=OUTDIR):
    alpha_distribution(df, identifiers, sns.stripplot, 'cluster', outdir=outdir)

def boxplot_alphas(df, identifiers, outdir=OUTDIR):
    alpha_distribution(df, identifiers, sns.boxplot, 'box', outdir=outdir)

def violin_alphas(df, identifiers, outdir=OUTDIR):
    alpha_distribution(df, identifiers, sns.violinplot, 'violin', outdir=outdir)

def swarm_alphas(df, identifiers, outdir=OUTDIR):
    alpha_distribution(df, identifiers, sns.swarmplot, 'swarm', outdir=outdir)

def category_dists_barplot(df, category, name, dist_col='neural_dist', outdir=OUTDIR):
    df = df[[category, dist_col]]
    cat_vals = []
    cat_means = []
    for cat_val, group in df.groupby(category):
        cat_vals.append(cat_val)
        cat_means.append(pylab.mean(group[dist_col]))
    order = pylab.argsort(cat_means)
    cat_vals = pylab.array(cat_vals)
    sorted_vals = cat_vals[order]
    pylab.figure()
    dist_plot = sns.barplot(x=category, y=dist_col, data=df, order=sorted_vals)
    pylab.xticks(rotation=90)#, size=20)
    if name == None:
        name = '%s/' % outdir
        if norm:
            name += 'norm_'
        name += 'pareto_dists_' + category + '.pdf'
    pylab.tight_layout()
    pylab.savefig(name + '_' + category + '.pdf', format='pdf')#, bbox_inches='tight')
    

def category_dists(df, categories, norm=False, outdir=OUTDIR):
    regression_df = df[['name', 'points', 'neural_dist']]
    regression_df = regression_df.drop_duplicates(subset='name')
    add_regression_cols(regression_df, 'points', 'neural_dist')
    regression_df['neural_dist_resid2'] = regression_df['neural_dist_resid'] ** 2
    for category in categories:
        df2 = df.drop_duplicates(subset=['name', category])
        df2 = remove_small_counts(df2, category)
        
        dist_col = ''
        if norm:
            dist_col += 'norm_'
        dist_col += 'neural_dist'
        
        name = '%s/' % OUTDIR
        if norm:
            name += 'norm_'
        name += 'pareto_dists'

        category_dists_barplot(df2, category, name, dist_col, outdir=outdir)

        df3 = regression_df.merge(df2, on='name')
        dist_col += '_resid2'
        name += '_regression'
        category_dists_barplot(df3, category, name, dist_col, outdir=outdir)


def scatter_dists(df, outdir=OUTDIR):
    df2 = df.drop_duplicates(subset='name')
    for neuron_type, group in df2.groupby('neuron_type'):
        print neuron_type, pylab.mean(group['neural_dist']), '+/-', pylab.std(group['neural_dist'], ddof=1)

    if LOG_DIST:
        df2 = df2[(df2['neural_dist'] > 0) & (df2['centroid_dist'] > 0) \
                 & (df2['random_dist'] > 0)]

#    df2 = df2[(df2['neural_dist'] > 1) & (df2['centroid_dist'] > 1) &\
#              (df2['random_dist'] > 1)]

    neural_dist = pylab.array(df2['neural_dist'])
    centroid_dist = pylab.array(df2['centroid_dist'])
    random_dist = pylab.array(df2['random_dist'])

    if LOG_DIST:
        neural_dist = pylab.log10(neural_dist)
        centroid_dist = pylab.log10(centroid_dist)
        random_dist = pylab.log10(random_dist)

    if CLIP_DIST:
        neural_dist = neural_dist.clip(min=0)
        centroid_dist = centroid_dist.clip(min=0)
        random_dist = random_dist.clip(min=0)

    assert len(neural_dist) == len(centroid_dist) == len(random_dist)

    order = pylab.argsort(neural_dist)
    neural_dist = neural_dist[order]
    centroid_dist = centroid_dist[order]
    random_dist = random_dist[order]

    x = range(len(neural_dist))
    pylab.figure()
    pylab.scatter(x, random_dist, c='m', label='random')
    pylab.scatter(x, centroid_dist, c='g', label='centroid')
    pylab.scatter(x, neural_dist, c='r', label='neural')
    pylab.xlim(-5, len(x))
    pylab.xlabel('neuron')
    ylab = ''
    if LOG_DIST:
        ylab += 'log-'
    ylab += 'distance'
    pylab.ylabel(ylab)
    #pylab.title('Distance to Pareto Front') 
    pylab.legend(loc='upper right')
    pylab.tight_layout()
    pylab.savefig('%s/pareto_dists.pdf' % outdir, format='pdf')#, bbox_inches='tight')
    pylab.close()

def alphas_hist(df, norm=False, outdir=OUTDIR):
    df2 = df.drop_duplicates(subset='name')
    alpha_col = ''
    if norm:
        alpha_col += 'norm_'
    alpha_col += 'alpha'
    alphas = list(df2[alpha_col])
    weights = pylab.ones_like(alphas) / len(alphas)
    pylab.figure()
    pylab.hist(alphas, range=(0, 1), weights=weights)
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    pylab.xlabel('alpha')
    pylab.ylabel('proportion')
    name = outdir + '/'
    if norm:
        name += 'norm_'
    name += 'alphas_hist.pdf'
    pylab.tight_layout()
    pylab.savefig(name, format='pdf')#, bbox_inches='tight')
    pylab.close()

def neuron_types_hist(df, norm=False, outdir=OUTDIR):
    df2 = df.drop_duplicates(subset='name')
    pylab.figure()
    alphas = []
    weights = []
    labels = []
    for neuron_type, group in df2.groupby('neuron_type'):
        alpha_col = ''
        if norm:
            alpha_col += 'norm_'
        alpha_col += 'alpha'
        alpha = list(group[alpha_col])
        #weight = pylab.ones_like(alphas) / len(alphas)
        weight = pylab.ones_like(alpha) / len(alpha)
        #pylab.hist(alphas, alpha=0.5, label=neuron_type,\
        #           range=(0, 1), weights=weights)
        alphas.append(alpha)
        weights.append(weight)
        labels.append(neuron_type)

    pylab.hist(alphas, range=(0, 1), weights=weights, label=labels)
    pylab.legend()
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    pylab.xlabel('alpha')
    pylab.ylabel('proportion')
    name = outdir + '/'
    if norm:
        name += 'norm_'
    name += 'neuron_types_hist.pdf'
    pylab.tight_layout()
    pylab.savefig(name, format='pdf')#, bbox_inches='tight')
    pylab.close()

def size_correlation(df):
    print "---------------size-alpha correlation----------------"
    df2 = df.drop_duplicates(subset='name')
    #df2['points'] = df2['points'].astype(int)
    corr1 = pearsonr(df2['points'], df2['alpha'])
    corr2 = spearmanr(df2['points'], df2['alpha'])
    print "pearson correlation: " + str(corr1) 
    print "spearman correlation: " + str(corr2)
    for name, group in df2.groupby('neuron_type'):
        print name
        corr1 = pearsonr(group['points'], group['alpha'])
        corr2 = spearmanr(group['points'], group['alpha'])
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
    print "random tree successes"
    print total_successes, total_trials
    print "p-value", float(total_successes) / total_trials

    df3 = df2[(df2['centroid_dist'] > 0) & (df2['random_dist'] > 0)]
    centroid_ratio = df3['neural_dist'] / df3['centroid_dist']
    random_ratio = df3['neural_dist'] / df3['random_dist']
    print "neural to centroid ratio", pylab.mean(centroid_ratio)
    print ttest_1samp(centroid_ratio, 1)
    print "neural to random ratio", pylab.mean(random_ratio)
    print ttest_1samp(random_ratio, 1)

    #print "dominate percentage", float(df2['dominates'].sum()) / df2['comparisons'].sum()
    #df3 = df2[df2['neuron_type'] != 'axon']
    #print "dendrite dominate percentage", float(df3['dominates'].sum()) / df3['comparisons'].sum()

    df4 = df2[df2['neural_dist'] < df2['centroid_dist']]
    centroid_trials = df2['neural_dist'].count()
    centroid_successes = df4['neural_dist'].count()
    assert centroid_trials >= centroid_successes
    print "beats centroid"
    print centroid_successes / float(centroid_trials)
    print binom_test(centroid_successes, centroid_trials)

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
        df2 = df2[df2['count'] >= 25]
        f = open(category + '.txt', 'w')
        with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
            print >> f,  df2
        f.close()

def neuron_type_alphas(df):
    df2 = df.drop_duplicates(subset=['name', 'neuron_type'])
    types = []
    dists = []
    indices = [0, 1, 2]
    for neuron_type, group in df2.groupby('neuron_type'):
        print neuron_type, pylab.mean(group['alpha']), '+/-', pylab.std(group['alpha'], ddof=1)
        types.append(neuron_type)
        dists.append(pylab.array(group['alpha']))
    for idx1, idx2 in combinations(indices, 2):
        type1, type2 = types[idx1], types[idx2]
        dist1, dist2 = dists[idx1], dists[idx2]
        print type1 + ' vs. ' + type2
        #print ttest_ind(dist1, dist2, equal_var=False)
        #print mannwhitneyu(dist1, dist2, alternative='two-sided')
        print ks_2samp(dist1, dist2)

def size_dist_correlation(df, outdir=OUTDIR):
    df2 = df.drop_duplicates(subset='name')
    df2 = df2[df2['neural_dist'] >= 1]
    points = df2['points']
    dists = df2['neural_dist']
    points = pylab.log10(points)
    dists = pylab.log10(dists) 
    print "---------size vs. dist correlation--------"
    print pearsonr(points, dists)
    print spearmanr(points, dists)

    add_regression_cols(df2, 'points', 'neural_dist', xtransform=pylab.log10,\
                        ytransform=pylab.log10)

    pylab.figure()
    pylab.scatter(points, dists)
    x = points
    y = df2['neural_dist_hat']
    order = pylab.argsort(x)
    x = x[order]
    y = y[order]
    pylab.plot(x, y, c='g')
    pylab.tight_layout()
    pylab.savefig('%s/size_dist.pdf' % outdir, format='pdf')#, bbox_inches='tight')
    pylab.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_file', default=OUTPUT_FILE)
    parser.add_argument('-c', '--categories_file', default=CATEGORIES_FILE)
    parser.add_argument('-m', '--models_file', default=MODELS_FILE)
    parser.add_argument('-f', '--figs_dir', default=OUTDIR)
    args = parser.parse_args()
    output_file = args.output_file
    categories_file = args.categories_file
    models_file = args.models_file
    assert data_file[-4:] == '.csv'
    outdir = args.outdir
    dfs = get_dfs()
    if TEST_NEW_FUNCTION:
        print dfs
        return None
    
    size_correlation(df)
    metadata(df)
    neuron_type_alphas(df)
    basic_stats(df)
    #categories_correlations(df)
    os.system('mkdir -p %s' % outdir)
    scatter_dists(df, outdir=outdir)
    boxplot_alphas(df, CATEGORIES, outdir=outdir)
    alphas_hist(df, norm=False, outdir=outdir)
    neuron_types_hist(df, norm=False, outdir=outdir)
    alphas_heat(df, CATEGORIES, outdir=outdir)
    dist_heats(df, CATEGORIES, DIST_FUNCS, outdir=outdir)
    category_dists(df, CATEGORIES, norm=False, outdir=outdir)
    size_dist_correlation(df, outdir=outdir)

    filtered_df = get_filtered_df(df)
    prefix = outdir[-4:]
    filtered_df.to_csv('%s_filtered.csv' % prefix, index=False)

if __name__ == '__main__':
    main()
