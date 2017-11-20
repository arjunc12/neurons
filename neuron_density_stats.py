import pareto_mst_plots
import neuron_density
import numpy as np
import pandas as pd
import pylab
from stats_utils import *

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
    pylab.figure()
    pylab.subplot(2, 1, 1)
    pylab.scatter(df['alpha'], np.log10(df['density']))
    pylab.subplot(2, 1, 2)
    pylab.scatter(df['alpha'], df['density_resid'])
    pylab.savefig('neuron_density/alpha_resid.pdf', format='pdf')
    pylab.close()

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
    pareto_mst_df = pareto_mst_plots.get_df()
    pareto_mst_df = pareto_mst_df[pareto_mst_df['cell_type'] != 'principal_cell']
    pareto_mst_df = pareto_mst_df[pareto_mst_df['neuron_type'] != 'axon']
    
    neuron_density_df = neuron_density.get_df()
    neuron_density_df.drop('points', axis=1, inplace=True)
    #neuron_density_df.drop_duplicates(subset='name', inplace=True)

    df = pd.merge(pareto_mst_df, neuron_density_df, on='name')
    df.drop_duplicates(subset='name', inplace=True)
    #print df[df['species'] == 'rabbit']
    df['density'] = df['mcost'] / df['volume']

    #print df['cell_type'].unique()
    #print df['species'].unique()
    #print df['region'].unique()

    add_regression_cols(df, 'mcost', 'density', np.log10, np.log10)
    for hue in [None, 'cell_type', 'species', 'region']:
        outliers = None
        if hue != None:
            outliers = biggest_outliers(df, hue)
            print hue, outliers
        make_density_plot(df, hue=hue, outliers=outliers)
    
    make_alpha_plots(df)

if __name__ == '__main__':
    main()
