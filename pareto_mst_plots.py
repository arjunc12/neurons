import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import pylab
import seaborn as sns
import os

OUTDIR = 'stats'

COLUMNS = ['name', 'cell_type', 'species', 'region', 'lab', 'alpha', 'neural_dist',\
            'centroid_dist', 'random_dist', 'trials', 'successes']

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
        dist_plot.tick_params(labelsize=5, axis='y')
        pylab.savefig('%s/%s_alphas_%s.pdf' % (OUTDIR, identifier, plot_descriptor), format='pdf')
        pylab.close

def cluster_alphas(df, identifiers):
    '''
    alphas = df['alpha']
    for identifier in identifiers:
        pylab.figure()
        stp = sns.stripplot(x='alpha', y=identifier, data=df, orient='h')
        #stp.set_xticklabels(stp.get_xticklabels(), rotation=75)
        stp.tick_params(labelsize=5)
        pylab.savefig('%s/%s_alphas.pdf' % (OUTDIR, identifier), format='pdf')
        pylab.close()
    '''
    alpha_distribution(df, identifiers, sns.stripplot, 'cluster')

def boxplot_alphas(df, identifiers):
    '''
    alphas = df['alpha']
    for identifier in identifiers:
        pylab.figure()
        bp = sns.boxplot(x='alpha', y=identifier, data=df, orient='h')
        #bp.set_xticklabels(bp.get_xticklabels(), rotation=75)
        bp.tick_params(labelsize=5)
        pylab.savefig('%s/%s_alphas_box.pdf' % (OUTDIR, identifier), format='pdf')
        pylab.close()
    '''
    alpha_distribution(df, identifiers, sns.boxplot, 'box')

def violin_alphas(df, identifiers):
    '''
    alphas = df['alpha']
    for identifier in identifiers:
        pylab.figure()
        vp = sns.violinplot(x='alpha', y=identifier, data=df, orient='h')
        #vp.set_xticklabels(bp.get_xticklabels(), rotation=75)
        vp.tick_params(labelsize=5)
        pylab.savefig('%s/%s_alphas_violin.pdf' % (OUTDIR, identifier), format='pdf')
        pylab.close()
    '''
    alpha_distribution(df, identifiers, sns.violinplot, 'violin')

def swarm_alphas(df, identifiers):
    '''
    alphas = df['alpha']
    for identifier in identifiers:
        pylab.figure()
        swp = sns.swarmplot(x='alpha', y=identifier, data=df, orient='h')
        #swp.set_xticklabels(bp.get_xticklabels(), rotation=75)
        swp.tick_params(labelsize=5)
        pylab.savefig('%s/%s_alphas_swarm.pdf' % (OUTDIR, identifier), format='pdf')
        pylab.close()
    '''
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
    pylab.scatter(x, pylab.log(neural_dist), c='r', label='neural')
    pylab.scatter(x, pylab.log(centroid_dist), c='b', label='centroid')
    pylab.scatter(x, pylab.log(random_dist), c='m', label='random')
    pylab.ylabel('log-distance')
    pylab.title('Distance to Pareto Front') 
    pylab.legend()
    pylab.savefig('%s/pareto_dists.pdf' % OUTDIR, format='pdf')
    pylab.close()

def main():
    fname = 'pareto_mst.csv'
    df = pd.read_csv(fname, names=COLUMNS)
    os.system('mkdir -p stats')
    #print df['species']
    scatter_dists(df)
    #cluster_alphas(df, ['species', 'cell_type', 'region'])
    boxplot_alphas(df, ['species', 'cell_type', 'region'])
    #violin_alphas(df, ['species', 'cell_type', 'region'])
    #swarm_alphas(df, ['species', 'cell_type', 'region'])

if __name__ == '__main__':
    main()
