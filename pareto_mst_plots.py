import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import pylab
import seaborn as sns
import os
from itertools import combinations

OUTDIR = 'stats'

'''
COLUMNS = ['name', 'cell_type', 'species', 'region', 'lab', 'alpha', 'neural_dist',\
            'centroid_dist', 'random_dist', 'trials', 'successes', 'comparisons',\
            'dominates']
'''
COLUMNS = ['name', 'cell_type', 'species', 'region', 'lab', 'points', 'alpha',\
           'neural_dist', 'centroid_dist', 'random_dist', 'trials', 'successes',\
            'comparisons', 'dominates']

NEURON_TYPE_LABELS = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite'}

def alphas_heat(df, categories):
    for cat1, cat2 in combinations(categories, 2):
        df2 = df.groupby([cat1, cat2], as_index=False).agg({'alpha' : pylab.mean})
        data = df2.pivot(cat1, cat2, 'alpha')
        pylab.figure()
        ax = sns.heatmap(data, vmin=0, vmax=1)
        pylab.savefig('%s/%s_%s_alphas_heat.pdf' % (OUTDIR, cat1, cat2), format='pdf')
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

def alpha_distribution(df, identifiers, plot_func, plot_descriptor):
    for identifier in identifiers:
        pylab.figure()
        dist_plot = plot_func(x='alpha', y=identifier, data=df, orient='h')
        dist_plot.tick_params(labelsize=10, axis='y')
        pylab.savefig('%s/%s_alphas_%s.pdf' % (OUTDIR, identifier, plot_descriptor), format='pdf')
        pylab.close()

def cluster_alphas(df, identifiers):
    alpha_distribution(df, identifiers, sns.stripplot, 'cluster')

def boxplot_alphas(df, identifiers):
    alpha_distribution(df, identifiers, sns.boxplot, 'box')

def violin_alphas(df, identifiers):
    alpha_distribution(df, identifiers, sns.violinplot, 'violin')

def swarm_alphas(df, identifiers):
    alpha_distribution(df, identifiers, sns.swarmplot, 'swarm')

def scatter_dists(df):
    neural_dist = df['neural_dist']
    centroid_dist = df['centroid_dist']
    random_dist = df['random_dist']
    
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
    pylab.ylabel('distance')
    pylab.title('Distance to Pareto Front') 
    pylab.legend()
    pylab.savefig('%s/pareto_dists.pdf' % OUTDIR, format='pdf')
    pylab.close()

def alphas_hist(df):
    alphas = list(df['alpha'])
    weights = pylab.ones_like(alphas) / len(alphas)
    pylab.hist(alphas, range=(0, 1), weights=weights)
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    pylab.savefig('%s/alphas_hist.pdf' % OUTDIR, format='pdf')

def neuron_types_hist(df):
    pylab.figure()
    for neuron_type, group in df.groupby('neuron_type'):
        alphas = list(group['alpha'])
        weights = pylab.ones_like(alphas) / len(alphas)
        pylab.hist(alphas, alpha=0.5, label=NEURON_TYPE_LABELS[neuron_type],\
                   range=(0, 1), weights=weights)
    pylab.legend()
    curr_ax = pylab.gca()
    curr_ax.set_ylim((0, 1))
    pylab.savefig('%s/neuron_types_hist.pdf' % OUTDIR, format='pdf')
    pylab.close()

def main():
    fname = 'pareto_mst.csv'
    df = pd.read_csv(fname, names=COLUMNS)
    df['neuron_type'] = df['name'].str[-1]
    df['neuron_type'] = df['neuron_type'].astype(int)
    #print df
    os.system('mkdir -p stats')
    #print df['species']
    scatter_dists(df)
    categories = ['species', 'cell_type', 'region']
    #cluster_alphas(df, categories)
    boxplot_alphas(df, categories)
    #violin_alphas(df, categories)
    #swarm_alphas(df, categories)
    alphas_hist(df)
    neuron_types_hist(df)
    alphas_heat(df, categories)

if __name__ == '__main__':
    main()
