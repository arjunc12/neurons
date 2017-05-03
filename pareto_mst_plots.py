import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import pylab
import seaborn as sns

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

def cluster_alphas(df, identifiers):
    alphas = df['alpha']
    for identifier in identifiers:
        sns.stripplot(x=identifier, y='alpha', data=df)
        pylab.savefig('figs/%s_alphas.pdf' % identifier, format='pdf')

def boxplot_alphas(df, identifiers):
    alphas = df['alpha']
    for identifier in identifiers:
        bp = sns.boxplot(x=df[identifier], y=df['alpha'])
        pylab.savefig('figs/%s_alphas_box.pdf' % identifier, format='pdf')

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
    pylab.savefig('figs/pareto_dists.pdf', format='pdf')
    pylab.close()

def main():
    columns = ['name', 'cell_type', 'species', 'region', 'lab', 'alpha', 'neural_dist',\
            'centroid_dist', 'random_dist', 'trials', 'successes']
    fname = 'pareto_mst.csv'
    df = pd.read_csv(fname, names=columns)
    scatter_dists(df)
    #cluster_alphas(df, ['species', 'cell_type', 'region'])
    #boxplot_alphas(df, ['species', 'cell_type', 'region'])
    

if __name__ == '__main__':
    main()
