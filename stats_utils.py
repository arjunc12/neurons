import numpy as np
import pandas as pd

def add_regression_cols(df, xcol, ycol, xtransform=None, ytransform=None):
    x = df[xcol]
    y = df[ycol]
    if xtransform != None:
        x = xtransform(x)
    if ytransform != None:
        y = ytransform(y)

    a, b = np.polyfit(x, y, 1)

    yhat = a * x + b
    yresid = y - yhat

    df[ycol + '_hat'] = yhat
    df[ycol + '_resid'] = yresid

