import matplotlib as mpl
mpl.use('agg')
import pylab
import pandas as pd
from pareto_steiner_stats import get_df, CATEGORIES, remove_small_counts

def pop_mean(df, column='alpha'):
    df2 = df.drop_duplicates(subset='name')
    return pylab.mean(df2[column])

def icc(df, category, mu=None):
    if mu == None:
        mu = pop_mean(df)

    df2 = df.drop_duplicates(subset=['name', category])
    betas = []
    errors = []
    for cat_val in df[category].unique():
        alphas = df2['alpha'][df2[category] == cat_val]
        alphas -= mu
        beta = pylab.mean(alphas)
        betas.append(beta)
        epsilon = alphas - beta
        errors += list(epsilon)

    beta_var = pylab.var(betas, ddof=1)
    error_var = pylab.var(errors, ddof=1)

    return beta_var / (beta_var + error_var)

def main():
    df = get_df()
    mu = pop_mean(df)
    for category in CATEGORIES:
        print category
        print icc(df, category, mu) 

if __name__ == '__main__':
    main()
