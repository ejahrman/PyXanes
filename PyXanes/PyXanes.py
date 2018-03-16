import os
import re
import itertools
import dill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from . import config
from .energyshift import energy_shift_angle


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
    '''warning, you're not allowed to save same name samples in different directories'''
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
    x, y = spectra.copy()
    return np.array([x, y/np.sum(y)])


class XesData:
    
    def __init__(self, sample, runnumber=None, batch=1, plotonload=True):
        '''Loads XES data taken with the RatTrap instrument.

        Args:
            sample: name of the sample as specified in the LabView
            runnumber: TODO
            batch: batch number specified in LabView

        TODO: do this better
        '''
        self.sample = sample
        self.runs = get_run_data(self.sample, dosum=False, runnumber=runnumber, batch=batch)
        self.spectrum = get_run_data(self.sample, dosum=True, runnumber=runnumber, batch=batch)
        if plotonload:
            quick_plot(sample,dosum=False,runnumber=runnumber,batch=batch)
            
    def plot(self, normalize=False, show=True, ax=None, shifted=False):
        # handle matplotlib ax for combined plots
        if ax is None:
            ax = plt.gca()

        # optional shift
        if shifted:
            spectrum = self.shiftedspectrum
        else:
            spectrum = self.spectrum
        
        #normalize, only 'integral' for now
        if normalize:#todo: peak normalization needed?
            spectrum = integral_normalize(spectrum)

        ax.plot(*spectrum, label=self.sample)

        if show:
            plotly_show()
            
    def save_to_file(self, filename):
        '''Save compressed version of this python object to the current directory'''
        print('Saving as: ', os.getcwd()+'\\'+filename)
        with open(filename,'wb') as f:
            dill.dump(self,f)
    
    def shift_spectrum(self, shiftenergy, eshift, yshift=0, yscale=1, normalize=False, resetanalyzer=False):
        if not hasattr(self, 'analyzerconfig') or resetanalyzer:
            print('Analyzer config not set, input now:')
            material = input('Analyzer material (si or ge)? ')
            h = int(input('[hkl] h? '))
            k = int(input('[hkl] k? '))
            l = int(input('[hkl] l? '))
            self.analyzerconfig = {
                'material': material,
                'h': h,
                'k': k,
                'l': l
            }
        self.shiftparams=dict(shiftenergy=shiftenergy,eshift=eshift,yshift=yshift,yscale=yscale,normalize=normalize)

        self.shiftedspectrum = energy_shift_angle(self.spectrum, shiftenergy, eshift, **self.analyzerconfig)

        if yshift != 0 or yscale != 1 or normalize:
            spectrum = self.shiftedspectrum
            if normalize:
                spectrum = integral_normalize(spectrum)
            x, y = spectrum
            self.shiftedspectrum = np.array([x, (y + yshift) * yscale])
        return self.shiftedspectrum

    def _add_lines_from_axes(self, ax=None):
        if ax is None:
            ax = plt.gca()
        datas = [l.get_data() for l in ax.lines]
        labels = [l.get_label() for l in ax.lines]
        plt.clf()
        return datas, labels

    def interact_shift(self, shiftenergy=None):
        from ipywidgets import interact, interactive
        import ipywidgets as widgets
        if shiftenergy is None:
            print('Unspecified shift energy, will use left-most datapoint.\nShifting from {} eV'.format(self.spectrum[0,0]))
            shiftenergy = self.spectrum[0,0]
        datas, labels = self._add_lines_from_axes()
        def compared_shifted(eshift, yshift, yscale, normalize=False):
            self.shift_spectrum(shiftenergy, eshift, yshift, yscale, normalize)
            self.plot(shifted=True, show=False)
            [plt.plot(*d, label=l) for (d,l) in zip(datas,labels)]
            plt.show()
        interact(
            compared_shifted, 
            eshift=widgets.FloatSlider(min=-10, max=10, step=0.01, value=0, continuous_update=True),
            yshift=widgets.FloatSlider(min=-0.01, max=0.01, step=0.0001, value=0, continuous_update=True),
            yscale=widgets.FloatSlider(min=0, max=2, step=0.001, value=1, continuous_update=True)
        )

    def difference(self):
        xf, yf = self.runs[max(self.runs.keys())]
        for k, (x, y) in self.runs.items():
            ydiff = yf - y
            if all(ydiff==0):
                plt.plot(x, ydiff, 'k')
            else:
                plt.plot(x, ydiff)
        plotly_show()

def load_XesData(filename):
    print('Attempting to load: ', os.getcwd()+'\\'+filename)
    with open(filename,'rb') as f:
        return dill.load(f)