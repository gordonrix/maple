#
#  DESCRIPTION   : Script for maple pipeline. Declares functions and classes that are used
#                   in several scripts throughout the pipeline
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import pathlib
import holoviews as hv

from bokeh.io import export_svgs
from selenium import webdriver as wd
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

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

def export_svg_plots(plots, file_name):
    """exports individual bokeh plots from a list of holoviews plots
    with the provided file name and the index of the plot in the list"""

    options = wd.ChromeOptions()
    options.add_argument('--headless')
    options.add_argument('--disable-gpu')
    options.add_argument("--no-sandbox")
    options.add_argument("--window-size=2000x2000")
    options.add_argument('--disable-dev-shm-usage')

    service = Service(ChromeDriverManager().install())
    webdriver = wd.Chrome(service=service, options=options)

    file_name_base = 'SVG_'+file_name[:-5]

    pathlib.Path(file_name_base).parent.absolute().mkdir(parents=True, exist_ok=True)
    
    for i, plot in enumerate(plots):
        fName = f'{file_name_base}_{i}.svg'
        p = hv.render(plot,backend='bokeh')
        p.output_backend='svg'
        export_svgs(p, 
            filename=fName,
            webdriver=webdriver)