import matplotlib as mpl
mpl.use('agg')
import pylab
import pandas as pd

OUTPUT_DIR = '/iblsn/data/Arjun/neurons/neuron_builder/snider/output'

PARAMETERS = {'snider': ['puncta radius', 'remove radius', 'trial length']}

def get_dfs():
    steiner_df = pd.read_csv('%s/pareto_steiner.csv' % OUTPUT_DIR, skipinitialspace=True)
    models_df = pd.read_csv('%s/models.csv' % OUTPUT_DIR, skipinitialspace=True)
    parameters_df = pd.read_csv('%s/parameters.csv' % OUTPUT_DIR, skipinitialspace=True)
    parameters_df = parameters_df[['tree'] + PARAMETERS['snider']]

    return steiner_df, models_df, parameters_df

def radius_plot(df):
    pylab.figure()
    pylab.scatter(df['radius'], df['alpha'])
    pylab.xlabel('radius')
    pylab.ylabel('alpha')
    pylab.savefig('steiner_stats/radius_alphas.pdf', format='pdf')
    pylab.close()

def alphas_distribution(df):
    pylab.figure()
    pylab.hist(df['alpha'])
    pylab.savefig('steiner_stats/neuron_builder_alphas.pdf', format='pdf')
    pylab.close()

def neuron_builder_stats(steiner_df, models_df):
    print pylab.mean(steiner_df['alpha'])
    df = models_df[models_df['model'] == 'neural']
    print pylab.mean(df['dist'])

def main():
    steiner_df, models_df, parameters_df = get_dfs()
    neuron_builder_stats(steiner_df, models_df)

if __name__ == '__main__':
    main()
