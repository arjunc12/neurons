import pareto_steiner_stats
import neuron_density
import numpy as np
import pandas as pd
import pylab
from stats_utils import *
from scipy.stats import pearsonr, spearmanr

def biggest_outliers(df, category, ycol='density'):
    resid_col = ycol + '_resid'
    max_val = None
    max_resid = 0
    for name, group in df.groupby(category):
        '''
        We are looking for a group that is consistently above/below the
        regression line. Thus we sum the residuals (instead of squared)
        residuals so that they will cancel out. We then take absolute value
        to measure magnitude.
        '''
        if len(group[resid_col]) < 5:
            continue
        resid = pylab.mean(group[resid_col])
        resid = abs(resid)
        if resid > max_resid:
            max_resid = resid
            max_val = [name]
        elif resid == max_resid:
            max_val.append(name)

    return max_val

def cat_to_color(cat_values):
    unique_values = set()
    cat_to_color = {}
    colors = []
    for value in cat_values:
        if value in unique_values:
            assert value in cat_to_color
        else:
            unique_values.add(value)
            cat_to_color[value] = len(unique_values)

        colors.append(float(cat_to_color[value]))

    colors = pylab.array(colors)
    #colors /= max(colors)
    return colors

def make_alpha_plots(df):
    alpha = df['alpha']
    density = df['density']
    log_density = np.log10(density)
    density_resid = df['density_resid']

    pylab.figure()
    neuron_types = df['neuron type'].unique()
    nrows = 3
    ncols = len(neuron_types) + 1
    enum_items = zip([density, log_density, density_resid],\
                     ['density', 'log-density', 'density residual'])
    for i, (y, label) in enumerate(enum_items):
        print "all neurons"
        #print label, pearsonr(alpha, y)
        print label, spearmanr(alpha, y)
        
        index1 = (i * ncols) + 1
        pylab.subplot(nrows, ncols, index1)
        pylab.scatter(alpha, y)
        pylab.ylabel(label)
        if i == 0:
            pylab.gca().set_title('all neurons')
        if i == len(enum_items) - 1:
            pylab.xlabel('alpha')

        for j, neuron_type in enumerate(neuron_types):
            x2 = alpha[df['neuron type'] == neuron_type]
            y2 = y[df['neuron type'] == neuron_type]
            print neuron_type
            #print label, pearsonr(x2, y2)
            print label, spearmanr(x2, y2)
            index2 = index1 + j + 1
            pylab.subplot(nrows, ncols, index2)
            pylab.scatter(x2, y2)
            if i == 0:
                pylab.gca().set_title(neuron_type)
            if i == len(enum_items) - 1:
                pylab.xlabel('alpha')
    
    pylab.tight_layout()

    pylab.savefig('neuron_density/alpha_resid.pdf', format='pdf')
    pylab.close()

    #print pearsonr(alpha, density)
    #print pearsonr(alpha, log_density)
    #print pearsonr(alpha, density_resid)


def sorted_points(df, col1, col2):
    idx = pylab.argsort(df[col1])
    x = df[col1][idx]
    y = df[col2][idx]
    return x, y

def make_density_plot(df, hue=None, outliers=None):
    c = 'b'
    cmap = None
    if hue != None:
        hue_values = df[hue]
        c = cat_to_color(hue_values)
        cmap = 'nipy_spectral'
    
    fname = 'neuron_density'
    if hue != None:
        fname += '_' + hue
    fname += '.pdf'

    pylab.figure()
    pylab.scatter(np.log10(df['mcost']), np.log10(df['density']), c=c, cmap=cmap)
    x, y = sorted_points(df, 'mcost', 'density_hat')
    pylab.plot(np.log10(x), y, c='r')

    if outliers != None:
        for outlier in outliers:
            df2 = df[df[hue] == outlier].copy()
            #print df2
            add_regression_cols(df2, 'mcost', 'density', np.log10, np.log10)
            x, y = sorted_points(df2, 'mcost', 'density_hat')
            pylab.plot(np.log10(x), y, c='g')

    pylab.savefig('neuron_density/%s' % fname)
    pylab.close()

def main():
    models_df, categories_df = pareto_steiner_stats.get_dfs()

    neuron_density_df = neuron_density.get_df()

    df = pd.merge(categories_df, neuron_density_df, on=['neuron name', 'neuron type', 'points'])
    
    df2 = df.drop_duplicates(subset=['neuron name', 'neuron type'])

    add_regression_cols(df2, 'volume', 'density', np.log10, np.log10)
    
    make_alpha_plots(df2)

if __name__ == '__main__':
    main()
