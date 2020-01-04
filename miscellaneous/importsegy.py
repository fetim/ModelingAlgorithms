#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:59:43 2019

@author: felipe
"""
#%%
import segyio
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np

def readSEGY(filename):
    print('Loading data cube from',filename,'with:')

    # Read full data cube
    data = segyio.tools.cube(filename)

    # Put temporal axis first
    data = np.moveaxis(data, -1, 0)

    #Make data cube fast to acess
    data = np.ascontiguousarray(data,'float32')

    #Read meta data
    segyfile = segyio.open(filename, "r")
    print('  Crosslines: ', segyfile.xlines[0], ':', segyfile.xlines[-1])
    print('  Inlines:    ', segyfile.ilines[0], ':', segyfile.ilines[-1])
    print('  Timeslices: ', '1', ':', data.shape[0])

    #Make dict with cube-info
    data_info = {}
    data_info['crossline_start'] = segyfile.xlines[0]
    data_info['inline_start'] = segyfile.ilines[0]
    data_info['timeslice_start'] = 1 #Todo: read this from segy
    data_info['shape'] = data.shape
    #Read dt and other params needed to do create a new


    return data, data_info

def parse_trace_headers(segyfile, n_traces):
    '''
    Parse the segy file trace headers into a pandas dataframe.
    Column names are defined from segyio internal tracefield
    One row per trace
    '''
    # Get all header keys
    headers = segyio.tracefield.keys
    # Initialize dataframe with trace id as index and headers as columns
    df = pd.DataFrame(index=range(1, n_traces + 1),
                      columns=headers.keys())
    # Fill dataframe with all header values
    for k, v in headers.items():
        df[k] = segyfile.attributes(v)[:]
    return df

if __name__ =="__main__":    

    filename =  "vp_marmousi-ii.segy"
    data_inline_full = []
    bin_headers_inline_full = []
    trace_headers_inline_full = []
    j = 0

    with segyio.open(filename,strict=False) as f:        
        n_traces = f.tracecount
        sample_rate = segyio.tools.dt(f) / 1000
        n_samples = f.samples.size
        twt = f.samples
        # Get all data into memory (could cause on big files)
        data_inline_full.append(f.trace.raw[:])  
        # Load headers
        bin_headers_inline_full.append(f.bin)
        trace_headers_inline_full.append(parse_trace_headers(f, n_traces))
    
    #%% Plot linha sismica

    print(np.shape(trace_headers_inline_full[0][:,0]))
    # line = 0
    # section = data_inline_full[line].T

    # vm = np.percentile(section, 95)
    # print("The 99th percentile is {:.0f}; the max amplitude is {:.0f}".format(vm, section.max()))
    # #%matplotlib inline
    # pl.figure(figsize=(18,6))
    # pl.imshow(section, cmap="jet", vmin=-vm, vmax=vm, aspect='equal')
    # pl.colorbar()
    # pl.show()