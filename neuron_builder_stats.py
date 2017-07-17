import matplotlib as mpl
mpl.use('agg')
import pylab
import pandas as pd

def radius_plot(df):
    pylab.figure()
    pylab.scatter(df['radius'], df['alpha'])
    pylab.xlabel('radius')
    pylab.ylabel('alpha')
    pylab.savefig('stats/radius_alphas.pdf', format='pdf')
    pylab.close()

def alphas_distribution(df):
    pylab.figure()
    pylab.hist(df['alpha'])
    pylab.savefig('stats/neuron_builder_alphas.pdf', format='pdf')
    pylab.close()

def neuron_builder_stats(df):
    print "mean alpha", pylab.mean(df['alpha'])
    print "mean dist", pylab.mean(df['dist'])

def main():
    df = pd.read_csv('neuron_builder.csv', names=['radius', 'size', 'dist', 'alpha'])
    radius_plot(df)
    alphas_distribution(df)
    neuron_builder_stats(df)

if __name__ == '__main__':
    main()
