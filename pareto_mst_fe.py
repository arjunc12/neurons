import matplotlib as mpl
mpl.use('agg')
import pandas as pd
from pareto_mst_plots import get_df, CATEGORIES, remove_small_counts
import pylab

def main():
    df = get_df()
    for category in CATEGORIES:
        df = remove_small_counts(df, category)
        df.drop('count', axis=1, inplace=True)

    regression_df = pd.DataFrame()
    for category in CATEGORIES:
        dummies = pd.get_dummies(df[category])
        regression_df = pd.concat([regression_df, dummies], axis=1)

    regression_df['name'] = df['name']

    regression_df = regression_df.groupby('name').agg(lambda x : int(pylab.any(x)))
    print regression_df


if __name__ == '__main__':
    main()
