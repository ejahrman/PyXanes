import os
import re
import itertools
import dill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from . import config


import plotly.offline as py
import plotly.tools as tls
py.init_notebook_mode()

def plotly_show():
    #get fig and convert to plotly
    fig = plt.gcf()
    plotlyfig = tls.mpl_to_plotly(fig, resize=True)
    
    #fix dumb automatic formatting choices
    plotlyfig['layout']['xaxis1']['tickfont']['size']=14
    plotlyfig['layout']['xaxis1']['titlefont']['size']=16
    plotlyfig['layout']['yaxis1']['tickfont']['size']=14
    plotlyfig['layout']['yaxis1']['titlefont']['size']=16
    plotlyfig['layout']['showlegend'] = True
    
    # plot
    py.iplot(plotlyfig)

def add_sample_dir_to_cache(sample):
    filere = re.compile(r'{}_\d+__EDX_\d+.txt'.format(sample))

    for root, dirs, files in os.walk(config.ROOTDIR):
        for file in files:
            if filere.match(file):
                config.DIRCACHE[sample] = '/'.join(root.split('\\')[:-1])
                break

def get_sample_path(sample):
    if sample in config.DIRCACHE:
        return config.DIRCACHE[sample]
    else:
        print('Directory for sample {} not yet cached, searching...'.format(sample))
        add_sample_dir_to_cache(sample)
        return config.DIRCACHE[sample]

def count_batches(path):
    '''counts number of subdirectories in the root directory, and assumes that each represents a batch'''
    return len([os.path.join(path, o) for o in os.scandir(path) if os.path.isdir(os.path.join(path,o))])

def load_sample(sample, runnumber=None, batch=1):
    '''
    By default, starts with run 0, and loads all files in order.
    
    If runnumber is 'int', then the runs are loaded starting from there.
    
    If runnumber is an interator, load each run in the iterator.
    '''
    if runnumber is None:
        runnumber = itertools.count()
    elif type(runnumber) is int:
        runnumber = itertools.count(runnumber)
    data = {}
    for rn in runnumber:
        try:
            file = '{0}/{1}_{2}/{1}{2}_alldata_{3}.txt'.format(get_sample_path(sample),sample,batch,rn)
            data[rn] = pd.read_table(file,header=2)
        except FileNotFoundError:
            break
    return data

def get_run_data(data, dosum=True, runnumber=None, batch=1):
    if type(data) is str:
        data = load_sample(data, runnumber=runnumber, batch=batch)
    spectra = {}
    for rn, df in data.items():
        spectra[rn] = df[['Energy_(eV)','cnts_per_live']].T.values
    if dosum:
        xsum, ysum = spectra.pop(rn)
        for rn, sp in spectra.items():
            x, y = sp
            assert all(xsum == x), 'Error, some data have different energies than others.'
            ysum = ysum + y
        spectra = np.array([xsum, ysum])
    return spectra

def quick_plot(sample, dosum=True, runnumber=None, batch=1, show=True):
    spectra = get_run_data(sample, dosum=dosum, runnumber=runnumber, batch=batch)
    if dosum:
        plt.plot(*spectra, label=sample)
    else:
        [plt.plot(*d, label=sample+str(i)) for i, d in spectra.items()]
    if show:
        plotly_show()


def integral_normalize(spectra):
    x, y = spectra
    return np.array([x, y/np.sum(y)])


class XesData:
    
    def __init__(self, sample, runnumber=None, batch=1, plotonload=True):
        self.sample = sample
        self.runs = get_run_data(self.sample, dosum=False, runnumber=runnumber, batch=batch)
        self.spectrum = get_run_data(self.sample, dosum=True, runnumber=runnumber, batch=batch)
        if plotonload:
            quick_plot(sample,dosum=False,runnumber=runnumber,batch=batch)
            
    def plot(self, normalize=None, show=True):
        if normalize == 'integral':
            plt.plot(*integral_normalize(self.spectrum), label=self.sample)
        else:
            plt.plot(*self.spectrum, label=self.sample)
        if show:
            plotly_show()
            
    def save_to_file(self, filename):
        '''Save compressed version of this python object to the current directory'''
        print('Saving as: ', os.getcwd()+'\\'+filename)
        with open(filename,'wb') as f:
            dill.dump(self,f)

def load_XesData(filename):
    print('Attempting to load: ', os.getcwd()+'\\'+filename)
    with open(filename,'rb') as f:
        return dill.load(f)