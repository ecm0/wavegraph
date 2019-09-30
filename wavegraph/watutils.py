# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

import os
import copy
import numpy as np

# import wavegraph_env as wgenv
# wgenv.init_ROOT()

import ROOT
from ROOT import wavearray
from ROOT import WDM
from ROOT import WSeries

import wseries
import timeseries

CWB_DEFAULT_INU = 6
CWB_DEFAULT_PRECISION = 10

def convert_numpyarray_to_wavearray(array, sampling_freq=1.):
    """
    Convert a Numpy array to a wavearray.
    
    Inputs:
    -------
    array   [Numpy Array] -- input data
    sampling_freq [float]-- sampling frequency [Hz]

    Output:
    -------
    out [Wavearray] -- Wavearray object
    """
    out = wavearray('double')(array.size)
    out.rate(sampling_freq)
    tmp_buffer = np.ndarray(shape=out.size(),
                            buffer=out.data,
                            dtype='double').T
    tmp_buffer[:] = array
    
    return out

def convert_wavearray_to_numpyarray(wavearray):
    """
    Convert a wavearray to a Numpy array

    Input:
    ------
    wavearray [Wavearray] -- input data

    Output:
    -------
    [Numpy Array] -- Numpy Array object
    """
    return np.ndarray(shape=wavearray.size(),
                      buffer=wavearray.data,
                      dtype='double').T

def initialize_WDM_types(grid, wat_type=None):
    """
    Return the WDM types associated to a CoherentWaveBurstGrid object

    Inputs:
    ------
    grid [CoherentWaveBurstGrid object] -- cWB grid
    wat_type                 [str/None] -- wavelet type 
    """
    if wat_type=='alternative_window':
        return [WDM('double')(int(2**scale),
                              'I',
                              int(2**scale), \
                              CWB_DEFAULT_INU, CWB_DEFAULT_PRECISION) \
                for scale in grid.timescales_exp]
    
    elif wat_type==None:
        return [WDM('double')(int(2**scale), int(2**scale), \
                              CWB_DEFAULT_INU, CWB_DEFAULT_PRECISION) \
                for scale in grid.timescales_exp]
    
def wdm_transform(signal, wdm_types, plotmode=False):
    """
    Return the WDM transforms and associated TF maps of the input signal

    Inputs:
    -------
    signal [Numpy Array] -- signal to decompose on WDM basis
    wdm_types     [list] -- WSeries type list

    Outputs:
    -------
    wdms    [list] -- WSeries objects list
    tfmaps  [list] -- list of 2D arrays of different shapes containing coefficients of 
                      WDM transform.
    """
    signal_wavearray = convert_numpyarray_to_wavearray(signal)

    wdms = []
    tfmaps = []
    time_axes = []
    freq_axes = []
    
    for wdm_type in wdm_types:
        
        wdm = WSeries('double')()
        wseries.extend_WSeries(wdm)
        wdm.Forward(signal_wavearray, wdm_type)
        wdms.append(wdm)
        direct, dual = wdm.nparray()
        
        tfmaps.append(dual**2 + direct**2)
        time_axes.append(wdm.times())
        freq_axes.append(wdm.freqs())

    if plotmode:
        return wdms, tfmaps, time_axes, freq_axes
    else:
        return wdms, tfmaps
    
def wilson_basis_func(wdm_type, num_samples, time_index, freq_index, dual_flag):
    """
    Returns the Wilson basis function that correspond to a given time and 
    frequency position.
    
    Inputs:
    -------
    wdm_type [WSeries] -- Wseries object
    num_samples  [int] -- number of samples for the reconstructed wavelet
    time_index   [int] -- time index of the pixel/wavelet to be reconstructed
    freq_index   [int] -- frequency index of the pixel/wavlet to be reconstructed
    dual_flag   [bool] -- enable/disable wavelet recosntruction on dual basis.
    
    Output
    ------
    out [Numpy Array] -- reconstructed wavelet.
    """
    # Compute wavelet associated with freq_coord.
    # time_coord does not matter here. This calls returns
    # a vector with the wavelet centered.
    # XXX test if time_coord can effectively be set to zero XXX
    wavelet = wavearray('double')()
    wdm_type.getBaseWave(int((wdm_type.m_Layer+1) * time_index + freq_index), \
                         wavelet, dual_flag)
    
    # Create a vector with num_samples samples and
    # copy the wavelet centered on time_index
    
    # First compute the indices where to copy the wavelet
    # Wrap around (circular mirroring) if part of
    # the wavelet support exceed vector boundaries
    wavelet_time = wdm_type.m_Layer * time_index
    idx_min = int(wavelet_time-(wavelet.size()-1)/2)
    idx_max = int(wavelet_time+(wavelet.size()-1)/2+1)
    idx = np.take(range(num_samples), range(idx_min,idx_max), mode='wrap')
    
    # Create output vector
    out = np.zeros(num_samples)
    
    # Copy the wavelet into the segment
    out[idx] = convert_wavearray_to_numpyarray(wavelet)
    
    return out

def load_nrms(filename, label=None):
    """  
    Load nRMS data from a .root or file. 
    The .txt or .txt.gz file cases are not yet implemented.

    Inputs:
    ------
    filename   [str] -- name of the .root, .txt or .txt.gz file to read.
    label [str/None] -- label of the leaf to read (.root file case).
    
    Output:
    -------
    nrms [ROOT.WSeries("double")/Numpy Array] -- nrms read in .root file
    """
    if not os.path.isfile(filename):
        raise Exception('File {} not found'.format(filename))
    
    if os.path.basename(filename).endswith('.root'):
        return ROOT.TFile(filename).Get(label)
    elif os.path.basename(filename).endswith(('.txt', '.txt.gz')):
        raise NotImplementedError
        #return np.loadtxt(filename)
    else:
        raise Exception('Unsupported format for nRMS files')

def whitening(*args):
    """
    Whitening in the time-frequency domain. 

    ASCII part currently not working.
    
    Input:
    ------
    signal                    [Timeseries object]  -- input signal
    noiserms [ROOT.WSeries("double")/Numpy Array] -- noise power spectral density
    wdm_type                            [Wseries] -- type of the WDM transform used for whitening

    Output:
    ------
    whitened signal [Timeseries object] -- whitened signal
    """
    if len(args) != 3:
        raise Exception("whitening() requires three input"
                        "args, {} given".format(len(args)))
    
    # Compute WDM transform
    wdm, _ = wdm_transform(args[0].data, args[2])
    wdm = wdm[-1] # squeeze list of 1 element

    wseries.extend_WSeries(wdm)
    
    if isinstance(args[1], ROOT.WSeries("double")):

        # Check consistency between noiserms and wdm_type
        if not args[1].maxLayer() == wdm.maxLayer():
            raise ValueError("The number of frequency bins of " \
                             "noiserms (={}) does not match that of the " \
                             "WDM transform (={})".format(args[1].maxLayer(),wdm.maxLayer()))
        
        # Whitening in the time-frequency domain
        wdm.white(args[1], +1)

    if isinstance(args[1], np.ndarray):

        direct, _ = wdm.nparray()
        
        # Check consistency between noiserms and wdm_type
        if direct.shape[0] != len(args[1]):
            raise ValueError("The number of frequency bins of " \
                             "noiserms (={}) does not match that of the " \
                             "WDM transform (={})".format(len(args[1]),wdm.maxLayer()))
    
        # Whitening: divide by amplitude noise spectrum, copy
        # the result "in place" (in the same array)
        np.copyto(direct, np.divide(direct.T, args[1]).T)
        
    # Invert the WDM transform to get whitened signal
    wdm.Inverse()
    
    # Convert from WSeries to Numpy array
    # and squeeze superfluous dimensions
    return timeseries.Timeseries(copy.copy(np.squeeze(wdm.nparray(False))),
                         args[0].sampling_freq,
                         0.0,
                         args[0].metadata)
