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
from scipy.stats import entropy, binom_test, ttest_1samp, ttest_ind, ttest_rel, wasserstein_distance, mannwhitneyu, wilcoxon, chisquare
from numpy.linalg import norm
import numpy as np
from stats_utils import *
import argparse
import neuron_density

from neuron_utils import DENDRITE_RATE
from check_robustness import INTERESTING_CELL_TYPES, INTERESTING_TRANSMITTERS

FIGS_DIR = 'steiner_stats'

TEST_NEW_FUNCTION = True
PAPER_PLOTS = False

OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output'

OUTPUT_FNAME = 'pareto_steiner.csv'
OUTPUT_FILE = '%s/%s' % (OUTPUT_DIR, OUTPUT_FNAME)

MODELS_FNAME = 'models.csv'
MODELS_FILE = '%s/%s' % (OUTPUT_DIR, MODELS_FNAME)

TRADEOFFS_FNAME = 'tradeoff_ratio.csv'
TRADEOFFS_FILE = '%s/%s' % (OUTPUT_DIR, TRADEOFFS_FNAME)

CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES_FILE_FILTERED = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories_filtered.csv'

CATEGORIES = ['cell type', 'species', 'region', 'neuron type', 'lab']

METADATA_DIR = '/iblsn/data/Arjun/neurons/metadata'

CATEGORY_MIN_COUNTS = {'species': 35, 'cell type' : 125, 'region' : 400, 'neuron type' : 100, 'lab' : 100}
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

def count_unique_neurons(df):
    return len(set(zip(df['neuron name'], df['neuron type'])))

def get_dfs(output_file=OUTPUT_FILE, categories_file=CATEGORIES_FILE,\
            models_file=MODELS_FILE, tradeoffs_file=TRADEOFFS_FILE):
    output_df = pd.read_csv(output_file, skipinitialspace=True)
    points_df = output_df[['neuron name', 'neuron type', 'points']]
    #print count_unique_neurons(output_df)

    models_df = pd.read_csv(models_file, skipinitialspace=True)
    models_df = pd.merge(models_df, points_df)
    #print count_unique_neurons(models_df)
    
    neural_df = models_df[models_df['model'] == 'neural']
    neural_df = neural_df[['neuron name', 'neuron type', 'dist']]

    cat_df = pd.read_csv(categories_file, skipinitialspace=True)
    categories_df = pd.merge(output_df, cat_df, on='neuron name')
    categories_df = pd.merge(categories_df, neural_df, on=['neuron name', 'neuron type'])
    #print count_unique_neurons(categories_df)

    tradeoffs_df = pd.read_csv(tradeoffs_file, skipinitialspace=True)
    tradeoffs_df = pd.merge(tradeoffs_df, points_df)
    #print count_unique_neurons(categories_df)

    density_df = neuron_density.get_df()
    #print count_unique_neurons(categories_df)

    models_df.drop_duplicates(inplace=True)
    categories_df.drop_duplicates(inplace=True)
    tradeoffs_df.drop_duplicates(inplace=True)
    density_df.drop_duplicates(inplace=True)

    return models_df, categories_df, tradeoffs_df, density_df

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

def dist_pval(counts1, counts2, dist_func):
    pass

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

        distribution1 = normalize_distribution(counts1)
        distribution2 = normalize_distribution(counts2)
        
        #distance1 = dist_func(counts1, counts2)
        #distance2 = dist_func(counts2, counts1)
        distance1 = dist_func(distribution1, distribution2)
        distance2 = dist_func(distribution2, distribution1)
        
        cat1 += [val1, val2]
        cat2 += [val2, val1]
        dist_vals += [distance1, distance2]

    dist_frame = pd.DataFrame()
    dist_frame[category + '1'] = cat1
    dist_frame[category + '2'] = cat2
    dist_frame['distance'] = dist_vals

    return dist_frame

def make_dist_frame_wasserstein(df, category):
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type', category])
    df2 = remove_small_counts(df2, category,\
                              min_count=CATEGORY_MIN_COUNTS[category])
 
    cat1 = []
    cat2 = []
    dist_vals = []
    for val1, val2 in combinations(list(df2[category].unique()), 2):
        sample1 = df2['alpha'][df2[category] == val1]
        sample2 = df2['alpha'][df2[category] == val2]
        dist1 = wasserstein_distance(sample1, sample2)
        dist2 = wasserstein_distance(sample2, sample1)
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
                   JSD : 'jsd', total_variation_distance: 'tvd',\
                   'wasserstein' : 'wasserstein'}

def dist_heat(df, category, alphas=None, dist_func=pseudo_kld, outdir=FIGS_DIR):
    dist_frame = None
    if dist_func == 'wasserstein':
        dist_frame = make_dist_frame_wasserstein(df, category)
    else:
        dist_frame = make_dist_frame(df, category, alphas, dist_func)
    
    
    dist_frame = dist_frame.pivot(category + '1', category + '2', 'distance')
    
    pylab.figure()
    
    ax = sns.heatmap(dist_frame, vmin=0, vmax=1)
    ax.set_ylabel(category + ' 1', fontsize=20)
    ax.set_xlabel(category + ' 2', fontsize=20)
    
    pylab.xticks(rotation='vertical', fontsize=20)
    pylab.yticks(rotation='horizontal', fontsize=20)
    #pylab.tight_layout()
    pylab.savefig('%s/%s_heat_%s.pdf' % (outdir, DIST_FUNC_NAMES[dist_func],\
                                         category.replace(' ', '_')),
                   format='pdf', bbox_inches='tight')
    pylab.close()

def kld_heat(df, category, alphas=None):
    dist_heat(df, category, alphas=alphas, dist_func=pseudo_kld)

def hellinger_heat(df, category, alphas=None):
    dist_heat(df, category, alphas=alphas, dist_func=hellinger_distance)

def jsd_heat(df, category, alphas=None):
    dist_heat(df, category, alphas=alphas, dist_func=JSD)

def wasserstein_heat(df, category, outdir=FIGS_DIR):
    dist_heat(df, category, alphas=None, dist_func='wasserstein', outdir=outdir)

def dist_heats(df, categories, dist_funcs, alphas=None, outdir=FIGS_DIR):
    for category in categories:
        for dist_func in dist_funcs:
            dist_heat(df, category, alphas=alphas, dist_func=dist_func, outdir=outdir)
        wasserstein_heat(df, category, outdir=outdir)

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

def val_distribution(df, val, categories, plot_func, plot_descriptor,\
                     outdir=FIGS_DIR, fig_suffix=None, category_subset=None,\
                     **kwargs):
    for category in categories:
        subset_cols = ['neuron name', 'neuron type', 'alpha']
        if category != 'neuron type':
            subset_cols.append(category)
        df2 = df.drop_duplicates(subset=subset_cols)

        if category_subset != None:
            df2 = df2[df2[category].isin(category_subset)]
        else:
            df2 = remove_small_counts(df2, category,\
                                      min_count=CATEGORY_MIN_COUNTS[category])

        cat_vals = []
        order_vals = []
        order_val = kwargs['order_val']
        for name, group in df2.groupby(category):
            cat_vals.append(name)
            if order_val == None:
                order_val.append(pylab.median(group[val]))
            else:
                order_vals.append(pylab.mean(group['dist']))
            print name, val, pylab.mean(group[val]), "+/-", pylab.std(group[val], ddof=1)
        
        cat_vals = pylab.array(cat_vals)
        #mean = pylab.array(medians)
        order = pylab.argsort(order_vals)
        order = cat_vals[order]

        pylab.figure()
        sns.set()
        dist_plot = plot_func(x=category, y=val, data=df2, order=order)
        dist_plot.tick_params(axis='x', labelsize=20, rotation=75)
        dist_plot.tick_params(axis='y', labelsize=20)
        #pylab.xlabel(category, fontsize=20)
        dist_plot.xaxis.label.set_visible(False)
        pylab.ylabel(val, fontsize=20)

        #sns.violinplot(x=category, y=val, data=df2, order=order)

        pylab.tight_layout()

        fname = '%s_%ss_%s' % (category.replace(' ', '_'),\
                              val.replace(' ', '_'), plot_descriptor)
        if fig_suffix != None:
            fname += '_%s' % fig_suffix
        pylab.savefig('%s/%s.pdf' % (outdir, fname), format='pdf')
        pylab.close()

def cluster_alphas(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                   category_subset=None, **kwargs):
    val_distribution(df, 'alpha', identifiers, sns.stripplot, 'cluster',\
                     outdir, fig_suffix, category_subset, **kwargs)

def boxplot_alphas(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                   category_subset=None, **kwargs):
    val_distribution(df, 'alpha', identifiers, sns.boxplot, 'box', outdir,\
                     fig_suffix, category_subset, **kwargs)

def violin_alphas(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                  category_subset=None, **kwargs):
    val_distribution(df, 'alpha', identifiers, sns.violinplot, 'violin',\
                     outdir, fig_suffix, category_subset, **kwargs)

def swarm_alphas(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                 category_subset=None, **kwargs):
    val_distribution(df, 'alpha', identifiers, sns.swarmplot, 'swarm', outdir,\
                     fig_suffix, category_subset, **kwargs)

def cluster_tradeoffs(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                      category_subset=None, **kwargs):
    val_distribution(df, 'tradeoff ratio', identifiers, sns.stripplot,\
                     'cluster', outdir, fig_suffix, category_subset, **kwargs)

def boxplot_tradeoffs(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                      category_subset=None, **kwargs):
    val_distribution(df, 'tradeoff ratio', identifiers, sns.boxplot, 'box',\
                     outdir, fig_suffix, category_subset, **kwargs)

def violin_tradeoffs(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                     category_subset=None, **kwargs):
    val_distribution(df, 'tradeoff ratio', identifiers, sns.violinplot,\
                     'violin', outdir, fig_suffix, category_subset, **kwargs)

def swarm_tradeoffs(df, identifiers, outdir=FIGS_DIR, fig_suffix=None,\
                    category_subset=None, **kwargs):
    val_distribution(df, 'tradeoff ratio', identifiers, sns.swarmplot,\
                     'swarm', outdir, fig_suffix, category_subset, **kwargs)
    
def category_dists(df, categories, outdir=FIGS_DIR, fig_suffix=None,\
                   category_subset=None, log_plot=True, **kwargs):
    for category in categories:
        df2 = df.drop_duplicates(subset=list(set(['neuron name', 'neuron type', category])))
        
        if category_subset != None:
            df2 = df2[df2[category].isin(category_subset)]
        else:
            df2 = remove_small_counts(df2, category,\
                                      min_count=CATEGORY_MIN_COUNTS[category])


        if category_subset != None:
            df2 = df2[df2[category].isin(category_subset)]

        if log_plot:
            df2['dist'] = pylab.log10(df2['dist'])
        
        
        cat_vals = []
        cat_means = []
        
        for cat_val, group in df2.groupby(category):
            dist = pylab.array(group['dist']).copy()

            cat_vals.append(cat_val)
            cat_mean = pylab.mean(dist)
            cat_means.append(cat_mean)

            if log_plot:
                dist = 10 ** dist
            print cat_val, pylab.mean(dist), "+/-", pylab.std(dist, ddof=1)

        order = pylab.argsort(cat_means)
        cat_vals = pylab.array(cat_vals)
        cat_means = pylab.array(cat_means)
        sorted_vals = cat_vals[order]
        sorted_means = cat_means[order]
        
        pylab.figure()
        sns.set()
        
        dist_plot = sns.barplot(x=category, y='dist', data=df2, order=sorted_vals)
        
        pylab.xticks(rotation=75, size=20)
        pylab.yticks(size=20)

        #pylab.xlabel(category, size=20)
        dist_plot.xaxis.label.set_visible(False)
        ylab = 'Distance to Pareto front'
        if log_plot:
            ylab = 'log(' + ylab + ')'
        pylab.ylabel(ylab, size=20)

        if 'ymin' in kwargs:
            pylab.ylim(ymin=kwargs['ymin'])
        if 'ymax' in kwargs:
            pylab.ylim(ymax=kwargs['ymax'])
        
        pylab.tight_layout()
        
        fname = 'pareto_dists_%s' % category.replace(' ', '_')
        if fig_suffix != None:
            fname += '_%s' % fig_suffix
        pylab.savefig('%s/%s.pdf' % (outdir, fname), bbox_inches='tight')

def scatter_dists(models_df, outdir=FIGS_DIR, subset=False):
    df = models_df[['neuron name', 'neuron type', 'model', 'dist']]
    model_dists = defaultdict(list)
    unique_neurons = 0
    for name, group in df.groupby(['neuron name', 'neuron type']):
        if len(group['model'].unique()) != 4:
            print group
            continue
        unique_neurons += 1
        if subset and unique_neurons % 5 != 0:
            continue
        for model, group2 in group.groupby('model'):
            #group2 = group2.head(n=20)
            model_dists[model].append(pylab.mean(group2['dist']))

    print "-------------"
    print "unique scatter neurons", unique_neurons
    print "-------------"
    
    order = pylab.argsort(model_dists['neural'])
    #x = pylab.arange(len(model_dists['neural']))
    pylab.figure()
    sns.set()
    
    model_colors = {'neural' : 'r', 'centroid' : 'g', 'random' : 'm', 'barabasi' : 'c'}
    model_markers = {'neural' : 'x', 'centroid' : 'o', 'random' : '^', 'barabasi' : 's'}
    model_labels = {'neural': 'Neural arbor', 'centroid' : 'Centroid', 'random' : 'Random', 'barabasi' : 'Barabasi-Albert'}
   
    max_dist = float('-inf')

    plot_order = ['random', 'barabasi', 'centroid', 'neural']
    for model in plot_order:
        dists = model_dists[model]
        dists = pylab.array(dists)
        y = dists[order]
        if LOG_DIST:
            y = pylab.log10(y)

        #y = y[::5]
        x = pylab.arange(len(y))

        max_dist = max(max_dist, max(y))

        color = model_colors[model]
        marker = model_markers[model]
        label = model_labels[model]
        pylab.scatter(x, y, label=label, c=color, marker=marker)
    
    pylab.xlabel('neural arbor index', fontsize=20)
    ylab = 'distance to Pareto front'
    if LOG_DIST:
        ylab = 'log(' + ylab + ')'
    pylab.ylabel(ylab, fontsize=20)
    leg = pylab.legend(ncol=2, frameon=True)
    leg.get_frame().set_linewidth(5)
    leg.get_frame().set_edgecolor('k')
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    pylab.ylim(0, max_dist + 0.6)

    pylab.xticks(fontsize=15)
    pylab.yticks(fontsize=15)

    fname = 'pareto_dists'
    if subset:
        fname += '_subset'
    
    pylab.tight_layout()
    pylab.savefig('%s/%s.pdf' % (outdir, fname), format='pdf')

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
        print "all neurons mean alpha", pylab.mean(alphas), "+/-", pylab.std(alphas, ddof=1)
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
    pylab.xlabel('alpha', size=20)
    pylab.ylabel('proportion', size=20)
    pylab.xticks(fontsize=20)
    pylab.yticks(fontsize=20)
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
        trials = len(group['success'])
        successes = sum(group['success'])
        #success_rate = pylab.mean(group['success'])

        ratios = group['ratio']
    
        '''
        successes = 0
        trials = 0
        ratios = []
        for (neuron_name, neuron_type), group2 in group.groupby(['neuron name', 'neuron type']):
            group2 = group2.head(n=20)
            successes += sum(group2['success'])
            trials += len(group2['success'])
            ratios.append(pylab.mean(group2['ratio']))
        '''

        print "success rate", float(successes) / float(trials), "trials", trials
        print "binomial p-value", binom_test(successes, trials)
        print "neural to %s ratio" % model, pylab.mean(ratios), "+/-", pylab.std(ratios, ddof=1)
        print "t-test p-value", ttest_1samp(ratios, popmean=1)

def null_models_check(models_df):
    df2 = models_df[models_df['model'] != 'neural']
    for (neuron_name, neuron_type), group in df2.groupby(['neuron name', 'neuron type']):
        ratios = group['ratio']
        models = group['model']
        r1 = list(ratios[models == 'centroid'])[:20]
        r2 = list(ratios[models == 'barabasi'])[:20]
        
        if pylab.mean(r1) > pylab.mean(r2):
            print neuron_name, neuron_type

def infmean(arr):
    return pylab.mean(masked_invalid(arr))

def metadata(df):
    print "-----------------------------------------------------"
    print "unique neurons"
    print len(set(zip(df['neuron name'], df['neuron type'])))

    for category in CATEGORIES:
        print "unique " + category
        print len(df[category].unique())
        df2 = df.drop_duplicates(subset=['neuron name', category])
        df2 = add_count_col(df2, category)
        #df2 = df2[df2['count'] >= 25]
        category_str = category.replace(' ', '_')
        f = open('%s/%s.txt' % (METADATA_DIR, category_str), 'w')
        with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
            print >> f,  df2
        f.close()

def neuron_type_alphas(df):
    print "-----------------------------------------------------"
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type'])
    types = []
    alphas = []
    dists = []
    for neuron_type, group in df2.groupby('neuron type'):
        print "------------"
        print neuron_type
        print "mean alpha", pylab.mean(group['alpha']), '+/-',\
                            pylab.std(group['alpha'], ddof=1)
        print "mean distance", pylab.mean(group['dist']), '+/-',\
                               pylab.std(group['dist'], ddof=1)
        types.append(neuron_type)
        alphas.append(pylab.array(group['alpha']))
        dists.append(pylab.array(group['dist']))
    indices = range(len(types))

    def alphas_chisquare(alphas1, alphas2):
        counts1 = defaultdict(int)
        counts2 = defaultdict(int)
        unique_alphas = set()
        for alpha in alphas1:
            counts1[alpha] += 1
            unique_alphas.add(alpha)
        for alpha in alphas2:
            counts2[alpha] += 1
            unique_alphas.add(alpha)

        f_obs = []
        f_exp = []
        for alpha in sorted(unique_alphas):
            f_obs.append(counts1[alpha])
            f_exp.append(counts2[alpha])

        return chisquare(f_obs, f_exp)

    for idx1, idx2 in combinations(indices, 2):
        print "------------"
        type1, type2 = types[idx1], types[idx2]
        alphas1, alphas2 = alphas[idx1], alphas[idx2]
        dists1, dists2 = dists[idx1], dists[idx2]
        print type1 + ' vs. ' + type2
        #print ttest_ind(dist1, dist2, equal_var=False)
        #print mannwhitneyu(dist1, dist2, alternative='two-sided')

        print "alphas ks-test", ks_2samp(alphas1, alphas2)
        print "alphas mann-whitney test", mannwhitneyu(alphas1, alphas2)
        print "alphas earth movers distance", wasserstein_distance(alphas1, alphas2)
        print "alphas chi square", alphas_chisquare(alphas1, alphas2)
        print "alphas welchs t-test", ttest_ind(alphas1, alphas2, equal_var=False)
        print "dists welch's t-test", ttest_ind(dists1, dists2, equal_var=False)

def vals_correlation(df, val1, val2, **kwargs):
    print "-----------------------------------------------------"
    print "%s-%s correlation" % (val1, val2)
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type'])
    
    if 'logtransform' in kwargs and kwargs['logtransform']:
        xtransform = pylab.log10
        ytransform = pylab.log10
        df2 = df2[(df2[val1] > 0) & (df2[val2] > 0)]
    
    v1 = df2[val1]
    v2 = df2[val2]
   
    xtransform = None
    ytransform = None

    if 'xtransform' in kwargs:
        xtransform = kwargs['xtransform']
    if 'ytransform' in kwargs:
        ytransform = kwargs['ytransform']
        
    if xtransform != None:
        v1 = xtransform(v1)
    if ytransform != None:
        v2 = ytransform(v2)

    print pearsonr(v1, v2)
    print spearmanr(v1, v2)

    regression_df = df2.copy()
    add_regression_cols(regression_df, val1, val2, xtransform=xtransform,\
                        ytransform=ytransform)

    grouping = None
    if 'grouping' in kwargs:
        grouping = kwargs['grouping']
    else:
        grouping = 'neuron type'
    
    grouping_subset = None
    if 'grouping_subset' in kwargs:
        assert 'grouping' in kwargs
        grouping_subset = kwargs['grouping_subset']
    else:
        #grouping_subset = ['axon', 'truncated axon', 'apical dendrite', 'basal dendrite']
        grouping_subset = df2[grouping].unique()
    
    sns.set()
    pylab.figure()
    nrows = len(grouping_subset) + 1
    pylab.subplot(nrows, 1, 1)
    pylab.scatter(v1, v2)
    x = v1
    y = regression_df['%s_hat' % val2]
    order = pylab.argsort(x)
    x = x[order]
    y = y[order]
    pylab.plot(x, y, c='g')

    row = 2

    df3 = df2[df2[grouping].isin(grouping_subset)]
    for name, group in df3.groupby(grouping):
        print name
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

    outdir = None
    if 'outdir' in kwargs:
        outdir = kwargs['outdir']
    else:
        outdir = FIGS_DIR
    figname = '%s/%s_%s.pdf' % (outdir, val1, val2)
    pylab.savefig('%s/%s_%s.pdf' % (outdir, val1, val2), format='pdf')

    pylab.close()

def size_dist_correlation(df, **kwargs):
    vals_correlation(df, 'points', 'dist', **kwargs)

def alpha_dist_correlation(df, **kwargs):
    vals_correlation(df, 'alpha', 'dist', **kwargs)

def size_alpha_correlation(df, **kwargs):
    vals_correlation(df, 'alpha', 'points', **kwargs)

def truncation_hist(df, outdir=FIGS_DIR):
    df2 = df[df['neuron type'].isin(['axon', 'truncated axon'])]
    df2 = df2.drop_duplicates(['neuron name', 'neuron type'])

    type_alphas = defaultdict(list)
    for neuron_name, group in df2.groupby('neuron name'):
        if len(group['neuron type']) < 2:
            continue
        for neuron_type, group2 in group.groupby('neuron type'):
            type_alphas[neuron_type] += list(group2['alpha'])

    alphas = []
    weights = []
    labels = []
    for neuron_type in type_alphas:
        alpha = type_alphas[neuron_type]
        alphas.append(alpha)
        weights.append(pylab.ones_like(alpha) / float(len(alpha)))
        labels.append(neuron_type)

    pylab.figure()
    sns.set()
    
    pylab.hist(alphas, range=(0, 1), weights=weights, label=labels)
    leg = pylab.legend(frameon=True)
    pylab.setp(leg.get_texts(), fontsize=20)
    leg_frame = leg.get_frame()
    leg_frame.set_linewidth(5)
    leg_frame.set_edgecolor('k')
    
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    
    pylab.xlabel('alpha', size=30)
    pylab.ylabel('proportion', size=30)
    pylab.tight_layout()
    
    name = 'truncation_hist'
    outname = '%s/%s.pdf' % (outdir, name)
    outname = outname.replace(' ', '_')
    pylab.savefig('%s/%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    axons = pylab.array(type_alphas['axon'])
    truncated_axons = pylab.array(type_alphas['truncated axon'])
    differences = axons - truncated_axons
    print '-------------------------'
    print "Truncation test"
    print min(differences), pylab.median(differences), max(differences)
    print pylab.mean(differences), "+/-", pylab.std(differences, ddof=1)
    print wilcoxon(axons, truncated_axons)
    print ttest_rel(axons, truncated_axons)

def triplet_analysis(df, categories=CATEGORIES):
    df2 = df.drop_duplicates()
    for category in categories:
        groupby_cols = categories[:]
        groupby_cols.remove(category)

        fname = 'triplets_%s.csv' % category
        fname = fname.replace(' ', '_')
        with open(fname, 'w') as f:
            for name, group in df2.groupby(groupby_cols):
                unique_vals = group[category].unique()
                group_items = ['--------------------', ', '.join(name), '--------------------']
                write_items = []
                for val1, val2 in combinations(group[category].unique(), 2):
                    sample1 = group['alpha'][group[category] == val1]
                    sample2 = group['alpha'][group[category] == val2]
                    n1 = len(sample1)
                    n2 = len(sample2)
                    if n1 > 25 and n2 > 25:
                        dist = wasserstein_distance(sample1, sample2)
                        write_items.append((dist, n1, n2, val1, val2))
                        #write_items.append('%s (%d), %s (%d), %f' % (val1, n1, val2, n2, dist))
                if len(write_items) > 0:
                    f.write('-----------------------\n')
                    f.write(', '.join(name) + '\n')
                    f.write('-----------------------\n')
                    write_items = reversed(sorted(write_items))
                    for dist, n1, n2, val1, val2 in write_items:
                        f.write('%s (%d), %s (%d), %f\n' % (val1, n1, val2, n2, dist))

def paired_categories_test(df, categories=CATEGORIES):
    for category in categories:
        print '----------------'
        print category
        print '----------------'
        subset_cols = ['neuron name', 'neuron type']
        if category != 'neuron type':
            subset_cols.append(category)
        df2 = df.drop_duplicates(subset=subset_cols)
        df2 = df2[subset_cols + ['alpha']]
        
        groupby_cols = subset_cols[:]
        groupby_cols.remove(category)

        for val1, val2 in combinations(df2[category].unique(), 2):
            print val1, val2
            sample1, sample2 = [], []
            for name, group in df2.groupby(groupby_cols):
                g1 = group[group[category] == val1]
                g2 = group[group[category] == val2]


                alpha1 = g1['alpha']
                alpha2 = g2['alpha']
                
                if len(alpha1) > 0 and len(alpha2) > 0:
                    s1 = pylab.mean(alpha1)
                    s2 = pylab.mean(alpha2)

                    sample1.append(pylab.mean(alpha1))
                    sample2.append(pylab.mean(alpha2))

            assert len(sample1) == len(sample2)
            if len(sample1) > 0:
                sample1 = pylab.array(sample1)
                sample2 = pylab.array(sample2)
                print "%s vs. %s" % (val1, val2)
                differences = sample1 - sample2
                print len(differences), pylab.mean(differences), pylab.std(differences, ddof=1)
                print ttest_rel(sample1, sample2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-od', '--output_dir', default=OUTPUT_DIR)
    parser.add_argument('-o', '--output_fname', default=OUTPUT_FNAME)
    parser.add_argument('-m', '--models_fname', default=MODELS_FNAME)
    parser.add_argument('-t', '--tradeoffs_fname', default=TRADEOFFS_FNAME)
    parser.add_argument('-c', '--categories_file', default=CATEGORIES_FILE)
    parser.add_argument('-f', '--figs_dir', default=FIGS_DIR)
    parser.add_argument('--synthetic', action='store_true')
    parser.add_argument('-r', '--rate', default=DENDRITE_RATE)
    parser.add_argument('--triplet', action='store_true')
    args = parser.parse_args()
    
    output_dir = args.output_dir
    output_fname = args.output_fname
    models_fname = args.models_fname
    tradeoffs_fname = args.tradeoffs_fname
    categories_file = args.categories_file
    figs_dir = args.figs_dir
    synthetic = args.synthetic
    triplet = args.triplet
    rate = float(args.rate)

    if synthetic:
        output_dir += '%0.1f' % rate
        output_fname = output_fname.replace('.csv', '_synthetic%0.1f.csv' % rate)
        models_fname = models_fname.replace('.csv', '_synthetic%0.1f.csv' % rate)
        tradeoffs_fname = tradeoffs_fname.replace('.csv', '_synthetic%0.1f.csv' % rate)
        figs_dir += '_synthetic'
        if rate != DENDRITE_RATE:
            figs_dir += '%0.1f' % rate
    
    os.system('mkdir -p %s' % figs_dir)

    output_file = '%s/%s' % (output_dir, output_fname)
    models_file = '%s/%s' % (output_dir, models_fname)
    tradeoffs_file = '%s/%s' % (output_dir, tradeoffs_fname)
    
    models_df, categories_df, tradeoffs_df, density_df = get_dfs(output_file,\
                                                                 categories_file,\
                                                                 models_file,\
                                                                 tradeoffs_file)
    
    categories_df = categories_df[categories_df['points'] >= MIN_POINTS]
    models_df = models_df[models_df['points'] >= MIN_POINTS]
    
    metadata(categories_df) 
    

    if TEST_NEW_FUNCTION:
        boxplot_alphas(categories_df, ['neuron type'], outdir=figs_dir, order_val='dist')
        boxplot_alphas(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='main',\
                       category_subset=INTERESTING_CELL_TYPES, order_val='dist')
        boxplot_alphas(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='transmitters',\
                       category_subset=INTERESTING_TRANSMITTERS, order_val='dist')

        violin_alphas(categories_df, ['region'], outdir=figs_dir, order_val='dist')
        return None
        
    if PAPER_PLOTS:
        neuron_type_alphas(categories_df)
        null_models_analysis(models_df)
        alphas_hist(categories_df, outdir=figs_dir)
    
        scatter_dists(models_df, outdir=figs_dir, subset=False)
        
        category_dists(categories_df, ['neuron type'], outdir=figs_dir, ymin=0, ymax=0.2)
        category_dists(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='main',\
                       category_subset=INTERESTING_CELL_TYPES)
        category_dists(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='transmitters',\
                       category_subset=INTERESTING_TRANSMITTERS)
       
        boxplot_alphas(categories_df, ['neuron type'], outdir=figs_dir, order_val='dist')
        boxplot_alphas(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='main',\
                       category_subset=INTERESTING_CELL_TYPES, order_val='dist')
        boxplot_alphas(categories_df, ['cell type'], outdir=figs_dir,\
                       fig_suffix='transmitters',\
                       category_subset=INTERESTING_TRANSMITTERS, order_val='dist')

        violin_alphas(categories_df, ['region'], outdir=figs_dir, order_val='dist')
    
        alpha_dist_correlation(categories_df, outdir=figs_dir)
        alpha_dist_correlation(categories_df, outdir=figs_dir, grouping='cell type', grouping_subset=INTERESTING_CELL_TYPES)
        alpha_dist_correlation(categories_df, outdir=figs_dir, grouping='cell type', grouping_subset=INTERESTING_TRANSMITTERS)
        return None
    
    print "-----------------------------------------------------"
    print "mean neural dist", pylab.mean(models_df['dist'][models_df['model'] == 'neural'])
    
    
    os.system('mkdir -p %s' % figs_dir)
     
    null_models_analysis(models_df)
    truncation_hist(categories_df, outdir=figs_dir)
        
    scatter_dists(models_df, outdir=figs_dir)
    
    neuron_types_hist(categories_df, outdir=figs_dir)
    alphas_hist(categories_df, outdir=figs_dir)
    
    alphas_heat(categories_df, CATEGORIES, outdir=figs_dir)
    dist_heats(categories_df, CATEGORIES, DIST_FUNCS, outdir=figs_dir)
    
    category_dists(categories_df, CATEGORIES, outdir=figs_dir)
    
    boxplot_alphas(categories_df, CATEGORIES, outdir=figs_dir)
    
    violin_tradeoffs(categories_df, CATEGORIES, outdir=figs_dir)
    
    
    size_dist_correlation(categories_df, outdir=figs_dir) 
    size_alpha_correlation(categories_df, outdir=figs_dir)
    
    cats = CATEGORIES[:]
    cats.remove('lab')
    if triplet:
        triplet_analysis(categories_df, cats)
        
    paired_categories_test(categories_df, categories=['neuron type'])
    
if __name__ == '__main__':
    main()
