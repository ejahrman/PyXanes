import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

    
def import_file(file):
    return pd.read_table(file, header=2)

def get_spectrum_from_dataframe(df):
    return df[['Energy_(eV)','cnts_per_live']].T.values

def plot_from_txt(file, show=True, plotcolumn='energy'):
    # possibility of multiple files to plot
    if type(file) != list:
        filelist = [file]
    else:
        filelist = file
    
    #set energy or angle to plot
    if plotcolumn == 'energy':
        colindex = 0
        xlabel = 'Energy (eV)'
    elif plotcolumn == 'angle':
        colindex = 1
        xlabel = 'Theta (deg)'
        
    #regex to get run-name from filepath
    namere = re.compile(r'.*/(.*)_alldata_(\d+).txt')
    
    for f in filelist:
        try:
            i = -1
            while True:
                i += 1
                data = get_spectrum_from_dataframe(import_file(f+'_alldata_{}.txt'.format(i)))
                name = namere.search(f+'_alldata_{}.txt'.format(i)).group(1)+'-'+namere.search(f+'_alldata_{}.txt'.format(i)).group(2)
                plt.plot(*data, label=name)
        except FileNotFoundError:
            continue
    plt.xlabel(xlabel)
    plt.ylabel('Counts per live second')
    if show:
        plotly_show()
		
def import_all(filelist, returnavg=True, removelist=[], verbose=False):
    if type(filelist) != list:
        filelist = [filelist]
    datas = []
    for f in filelist:
        try:
            i = -1
            while True:
                i += 1
                if i in removelist:
                    continue
                data = get_spectrum_from_dataframe(import_file(f+'_alldata_{}.txt'.format(i)))
                if verbose:
                    print('Imported: ',f+'_alldata_{}.txt'.format(i))
                datas.append(data)
        except FileNotFoundError:
            continue
    datas = np.array(datas)
    if returnavg:
        x = datas[0,0]
        for d in datas:
            assert (d[0] == x).all(), 'Not all runs have same energy values'
        datas = np.array([x,np.mean(datas[:,1],axis=0)])
    return datas
	
def mu(izero, it):
    xz, yz = izero
    xt, yt = it
    assert (xz == xt).all(), 'Need to have same energy values'
    return np.array([xz, np.log(yz/yt)])
	
def check_scan_progress(path='./NiIt551Extended/NiIt551Extended_1/', root='NiIt551Extended_', index=2):
    numedxs = len(glob.glob(path+root+'{}*'.format(0)))
    numcurr = len(glob.glob(path+root+'{}*'.format(index)))
    print( '{} datapoints left to go!\r'.format(numedxs-numcurr), end='')
#     print( '|'+'-'*int((numcurr/numedxs*100))+'>|\r'.format(numedxs-numcurr), end='')