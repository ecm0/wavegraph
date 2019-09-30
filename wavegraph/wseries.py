# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

# Extension of the WSeries class that adds few utilities

import os

import numpy as np

# import wavegraph_env as wgenv
# wgenv.init_ROOT()

from ROOT import WSeries
from ROOT import wavearray

def extend_WSeries(wseries):
    """ Changes the class of wseries (a wseries object) 
    to WSeriesExtended. """
    wseries.__class__ = WSeriesExtended
    
class WSeriesExtended(WSeries("double")):
    """ Extend a WSeries object in order to extract
    the informations it contains. """
    def num_times(self):
        """ Return the number of time bins """
        return self.sizeZero()

    def num_freqs(self):
        """ Return the number of frequency bins """
        return self.maxLayer()+1
    
    def num_samples(self):
        """ Return the total number of samples """
        return self.pWavelet.nWWS

    def time_sampling(self):
        """ Return the time sampling step in sec """
        if self.resolution() > 0:
            return 1.0/(2 * self.resolution())
        else:
            return None

    def freq_sampling(self):
        """ Return the frequency sampling step in Hz """
        return self.resolution()
    
    def times(self):
        """ Return the time axis """
        return self.time_sampling() * np.arange(self.num_times())

    def freqs(self):
        """ Return the frequency axis """
        return self.freq_sampling() * np.arange(self.num_freqs())
    
    def nparray(self, with_dual=True):
        """
        Return the WSeries contents in a Numpy array.

        with_dual -- If true, returns a pair (direct array, dual
        array). Otherwise, only returns the direct array.
        """
        # This is how the wavelet data are stored
        #      DataType_t *pWWS;     //! pointer to wavelet work space      
        #      unsigned long nWWS;   // size of the wavelet work space

        if self.pWavelet.nWWS == 0:
            return np.empty(0)

        direct_array = np.ndarray(
            shape=(self.num_times(), self.num_freqs()),
            buffer=self.pWavelet.pWWS, dtype='double').T

        # These are other possible mappings from double* to NumPy arrays
        #     arr = numpy.frombuffer(buff,count=N)
        #     arr = numpy.array(buff,copy=True)

        if with_dual:
            dual_array = np.ndarray(
                shape=(self.num_times(),self.num_freqs()),
                buffer=self.pWavelet.pWWS,
                # XXX Why float_() instead of double??
                offset=(self.maxIndex()+1)*np.float_().itemsize,
                dtype='double').T

        return (direct_array, dual_array) if with_dual else direct_array
    
    def wavearray(self):
        """ Return the WSeries contents in a wavearray. """

        # This is how the wavelet data are stored
        #      DataType_t *pWWS;     //! pointer to wavelet work space      
        #      unsigned long nWWS;   // size of the wavelet work space

        if self.pWavelet.nWWS == 0:
            return wavearray('double')(0)
        
        output = wavearray('double')(self.pWavelet.pWWS, self.pWavelet.nWWS, 0)
        output.rate(self.rate())

        return output
