# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

from __future__ import division

import os
import sys
import h5py
import numpy
import argparse
import getpass
import time
import logging
import fractions
import subprocess
import sys
from scipy import signal

LOGGING_FMT = "%(levelname)s -- %(filename)s:line %(lineno)s in %(funcName)s(): %(message)s"


class Timeseries(object):
    def __init__(self, data, fs, t0=0, meta=None):
        self.data = numpy.array(data)
        self.sampling_freq = fs
        self.t0 = t0
        self.metadata = meta
        
    def duration(self):
        return len(self.data)/self.sampling_freq
    
    def time(self):
        return self.t0 + numpy.arange(len(self.data))/self.sampling_freq

    def __str__(self):
        return str(self.data)
    
def read_timeseries(filename):
    """
    Read an .hdf file and returns a list of Timeseries.
    
    Input
    -----
    filename -- input .hdf file [string]
    
    Output
    -----
    timeseries [list] -- list of Timeserie objects
    infos [str] -- description of the data set.
    
    .hdf5 file structure
    --------------------
    - file
    --- description [string]
    --- metadata [string]
    - data
    - timeseries #0 [Numpy array]
    --- sampling_freq [float, Hz]
    --- descr [str]
    - ...
    
    """
    if not os.path.isfile(filename):
        raise Exception('File {} not found'.format(filename))
    
    logging.info("Reading {}".format(filename))
        
    with h5py.File(filename, 'r') as file:
        timeseries = [Timeseries(numpy.array(ts), \
                                 ts.attrs['sampling_freq'], \
                                 0.0, \
                                 ts.attrs['descr']) \
                      for name, ts in file['/data'].items()]
        
        infos = {'description:': file.attrs['description'],
                 'metadata': file.attrs['metadata']}

    return timeseries, infos

def write_timeseries(filename, timeseries, metadata):
    """
    Write timeseries to a .hdf file.
    
    Input
    -----
    filename [str]  -- name of .hdf file 
    timeseries [list] -- iterable of Timeseries objects
    metadata [str] -- additional info about input data
    
    .hdf5 file structure
    --------------------
    - file
    --- description [string]
    --- metadata [string]
    - data
    - timeseries #0 [Numpy array]
    --- sampling_freq [float, Hz]
    --- descr [string]
    - ...
    """
    if os.path.isfile(filename):
        logging.info('{} already exists -- Removing'.format(filename))
        os.remove(filename)

    # Open writing session for .hdf file.
    try:
        outfile = h5py.File(filename, 'w')
    except IOError:
        raise IOError('Cannot write file {}'.format(filename))
        
    # Create header...
    try:
        githash = subprocess.check_output(["git", "describe", "--always"])
    except OSError:
        logging.warning("Git executable not found -- Unable to get Git hash")
        githash = "unknown"
        
    outfile.attrs['description'] = \
                'Generated by {} at {} -- Git: {}'.format(getpass.getuser(),
                                                          time.strftime('%Y-%m-%d %H:%M:%S'),
                                                          githash)
    outfile.attrs['metadata'] = metadata
    
    # ... then store attributes and timeseries.
    group = outfile.create_group('data')
    for n, ts in enumerate(timeseries):
        
        dataset = group.create_dataset('timeseries #{}'.format(n), data=ts.data, \
                                       shape=ts.data.shape, compression='gzip')
        dataset.attrs['sampling_freq'] = ts.sampling_freq
        dataset.attrs['descr'] = ts.metadata

    # Close writing session of .hdf file.
    outfile.close()
        
    logging.info('Wrote {} timeseries in {}'.format(len(timeseries), filename))

## XXX Could be included in the above class XXX
def resample(s, p, q, h=None):
    """Change sampling rate by rational factor. This implementation is based on
    the Octave implementation of the resample function. It designs the 
    anti-aliasing filter using the window approach applying a Kaiser window with
    the beta term calculated as specified by [2].
    
    Ref [1] J. G. Proakis and D. G. Manolakis,
    Digital Signal Processing: Principles, Algorithms, and Applications,
    4th ed., Prentice Hall, 2007. Chap. 6
    Ref [2] A. V. Oppenheim, R. W. Schafer and J. R. Buck, 
    Discrete-time signal processing, Signal processing series,
    Prentice-Hall, 1999
    """
    
    gcd = fractions.gcd(p,q)
    if gcd > 1:
        p = p/gcd
        q = q/gcd
        
    if h is None: #design filter
        
        #properties of the antialiasing filter
        log10_rejection = -3.0
        stopband_cutoff_f = 1.0/(2.0 * max(p,q))
        roll_off_width = stopband_cutoff_f / 10.0
        
        #determine filter length
        #use empirical formula from [2] Chap 7, Eq. (7.63) p 476
        rejection_db = -20.0*log10_rejection;
        l = numpy.ceil((rejection_db-8.0) / (28.714 * roll_off_width))
        
        #ideal sinc filter
        t = numpy.arange(-l, l + 1)
        ideal_filter=2 * p * stopband_cutoff_f * numpy.sinc(2*stopband_cutoff_f*t)
        
        #determine parameter of Kaiser window
        #use empirical formula from [2] Chap 7, Eq. (7.62) p 474
        beta = signal.kaiser_beta(rejection_db)
        
        #apodize ideal filter response
        h = numpy.kaiser(2*l+1, beta) * ideal_filter
        
    ls = len(s)
    lh = len(h)
    
    l = (lh - 1)/2.0
    ly = int(numpy.ceil(ls*p/float(q)))
    
    #pre and postpad filter response
    nz_pre = int(numpy.floor(q - numpy.mod(l,q)))
    hpad = h[-lh+nz_pre:]
    
    offset = int(numpy.floor((l+nz_pre)/q))
    nz_post = 0;
    
    while numpy.ceil(((ls-1)*p + nz_pre + lh + nz_post )/q ) - offset < ly:
        nz_post += 1
        
    hpad = hpad[:lh + nz_pre + nz_post]
    
    #filtering
    xfilt = upfirdn(s, hpad, p, q)
    
    return xfilt[offset-1:offset-1+ly], h

def upfirdn(s, h, p, q):
    """
    Upsample signal s by p, apply FIR filter as specified by h, and 
    downsample by q. Using fftconvolve as opposed to lfilter as it does not seem
    to do a full convolution operation (and its much faster than convolve).
    """

    a = upsample(s, p)
    b = signal.fftconvolve(h, a)
    return downsample(b, q)