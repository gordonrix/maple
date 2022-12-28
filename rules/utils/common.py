#
#  DESCRIPTION   : Script for maple pipeline. Declares functions and classes that are used
#                   in several scripts throughout the pipeline
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np

def dist_to_DF(dist, x, y):
    """given a np.bincount output, i.e. a distribution of values in which each value is
    the number of observations corresponding to the value's position in the distribution
    (e.g. [0,3,1] is a distribution for 4 values in which 3 of the values were '1'),
    this will calculate the distribution as a proportion and as a cumulative proportion
    and produce a DataFrame from these three numpy arrays that describe the distribution
    
    dist:       np.array of shape (maximumValue)
    x, y:       strings that describes the x and y variable"""

    maxVal = dist.shape[0]
    dist = dist.reshape(1,maxVal)

    values = np.arange(maxVal).reshape(1,maxVal)
    proportion = np.divide(dist, dist.sum())
    cumsum = np.cumsum(proportion).reshape(1,maxVal)

    df = pd.DataFrame( np.concatenate((values, dist, proportion, cumsum), axis=0).T, columns = [x, f"total {y}", f"proportion of {y}", f"cumulative proportion of {y}"] )
    df[df.columns[[0,1]]] = df[df.columns[[0,1]]].apply(lambda x: x.astype(int)) # convert x value and total counts columns to int

    return df