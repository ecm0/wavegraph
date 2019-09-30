# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

import numpy
import os
import itertools
import collections
import h5py
import getpass
import time
import subprocess
import logging

class CoherentWaveBurstGrid(object):
    """
    Time-frequency-scale grid associated to the WDM transform used by
    the coherent WaveBurst pipeline.

    Main attributes:
    ----------------
    sampling_freq -- sampling frequency; numbers of samples of the signal per second
    timescales -- a list of increasing timescales associated with each grid plane
    timescales_exp -- a list of increasing timescale exponents associated with each grid plane
    timescale_min, timescale_max -- minimum and maximum timescales
    """
    def __init__(self, sampling_freq, min_scale_exp, max_scale_exp):
        """
        Main attributes:
        ----------------
        sampling_freq -- sampling frequency; numbers of samples of the signal per second
        timescales -- a list of increasing timescales associated with each grid plane                                             
        timescales_exp -- a list of increasing timescale exponents associated with each grid plane
        timescale_min, timescale_max -- minimum and maximum timescales
        """
        self.sampling_freq = float(sampling_freq)
        self.timescales_exp = numpy.arange(min_scale_exp, max_scale_exp+1)
        self.timescales = 2**numpy.arange(min_scale_exp,
                                          max_scale_exp+1)/self.sampling_freq
        
        # Convenient shortcuts:
        self.timescale_max = self.timescales[-1]
        self.timescale_min = self.timescales[0]
        
        
class GridPoint(collections.namedtuple('GridPoint', 'scale_index time_index freq_index')):
    """ Point or pixel in a CoherentWaveBurstGrid. """
    
    def phys_coords(self, grid):
        return {"scale": grid.timescales_exp[self.scale_index], \
                "time": self.time_index * grid.timescales[self.scale_index], \
                "frequency": self.freq_index \
                / (2 *grid.timescales[self.scale_index])}

    @classmethod
    def from_phys_coords(cls, grid, timescale, time, freq):
        """
        Instantiate a GridPoint object from physical coordinates
        """
        scale_index = numpy.argwhere(grid.timescales == timescale)
    
        if not scale_index:
            logging.warning('Requested timescale does not exist in grid')
            return None

        return cls(scale_index, \
                         int(times / timescale), \
                         int(2 * freq * timescale))

ClusterNamedTuple = collections.namedtuple(
    "Cluster", "grid_points values")

class Cluster(ClusterNamedTuple):
    """
    Cluster on a CoherentWaveBurstGrid, with values associated to its
    nodes (the values are defined by the user).
    The instance attributes are those of ClusterNamedTuple.
    Cluster objects are hashable. This can be useful for
    avoiding duplicate clusters in the grid.
    """
    def __new__(cls, grid_points, values, metadata):
        """
        grid_points -- iterable of GridPoints.
        values -- iterable of associated values (of the same size as
        grid_points).
        metadata -- description string
        """
        instance = super(Cluster, cls).__new__(cls,
                                               # Hashable components:
                                               tuple(grid_points), tuple(values))
        instance.metadata = metadata
        return instance

    @classmethod
    def from_numpyarray(cls, array, metadata):
        """
        array -- Numpy named array with fields: 'scale_index', 'time_index',
                'freq_index' and 'value'
        """
        return cls(tuple(GridPoint(p['scale_index'], \
                                   p['time_index'],  \
                                   p['freq_index'])  \
                         for p in array),  \
                   tuple(array['value']), \
                   metadata)
    
    def to_numpyarray(self):
        """ Convert a Cluster object into a Numpy array object. """
        out = [(p.scale_index, p.time_index, p.freq_index, v) \
               for p, v in zip(self.grid_points, self.values)]
        return numpy.array(out,
                           dtype=[('scale_index', int), ('time_index', int), \
                                  ('freq_index', int), ('value', float)])
    
    def phys_coords(self, grid):
        """ 
        Convert list of GridPoint objects into a list of 
        physical quantities. 
        """
        phys_points = [p.phys_coords(grid) for p in self.grid_points]
        return numpy.array([(p['scale'], p['time'], p['frequency']) for p in phys_points],
                           dtype=[('scale', int), ('time', float), ('frequency', float)])

    def reject_pixels_at_zero_freq(self):
        """ Remove pixels whose frequency dimension is zero. """
        indices = self.to_numpyarray()['freq_index'] != 0
        return Cluster(self.grid_points[indices], self.values[indices])

    def phys_values(self, grid):
        """ 
        Convert values from cluster of indices to
        a cluster of physical quantities.
        """
        phys_values = [numpy.sqrt(val * grid.timescales[gp.scale_index]) \
                       for (val, gp) in zip(self.values, self.grid_points)]
        return numpy.array(phys_values)

def shift_clusters_to_zero_index(clusters, grid):
    """ 
    Shift all computed clusters to zero time index.

    clusters [Numpy Array] -- list of Cluster objects
    grid [CoherentWaveBurstGrid object] -- cWB grid
    """
    # Find the maximum occupied scale. This is needed because time
    # shifts on the grid can only be done in steps of the largest scale:
    scale_index_max = max(
        node.scale_index for cluster in clusters for node in cluster.grid_points)

    # Project the smallest time, at any scale, onto its time index
    # at the maximum occupied scale:  # XXX What for? algorithm?
    time_min = min(node.time_index*grid.timescales[node.scale_index]
                   for cluster in clusters for node in cluster.grid_points)

    index_shift_at_scale_max = int(numpy.floor(
        time_min/grid.timescales[scale_index_max]))
    
    return [
        Cluster(
            [grid_point._replace(time_index=(
                grid_point.time_index
                - index_shift_at_scale_max
                *2**(scale_index_max - grid_point.scale_index)))
             for grid_point in cluster.grid_points],
        
            cluster.values,
            cluster.metadata
        )
            for cluster in clusters]

def read_clusters(filename):
    """
    Read clusters from a .hdf file.

    Input
    -----
    filename -- name of .hdf file [str]

    .hdf5 file structure
    --------------------
    - file
    --- description [str]
    --- metadata [str]
    - data
    - cluster #0 [Numpy array]
    --- descr [str]
    - ...
 
    """
    if not os.path.isfile(filename):
        raise Exception('File {} not found'.format(filename))
    
    logging.info("Reading {}".format(filename))
    
    with h5py.File(filename, 'r') as file:
        clusters = [Cluster.from_numpyarray(numpy.array(c), \
                                    c.attrs['description']) \
                    for name, c in file['/data'].items()]
        infos = {'description': file.attrs['description'], \
                 'metadata': file.attrs['metadata'], \
                 'params': file.attrs['params']}

    return clusters, infos

def write_clusters(filename, clusters, params, metadata):
    """
    Write clusters to a .hdf file.

    Input
    -----
    filename -- name of .hdf file [str]
    clusters -- iterable of Cluster objects
    params   -- parsable string with main parameters [str]
    metadata -- additional info about the input data [str]

    .hdf5 file structure
    --------------------
    - file
    --- description [str]
    --- metadata [str]
    --- params [dict]
    - data
    - cluster #0 [Numpy array]
    --- descr [str]
    - ...
    """
    # Open writing session for .hdf file.
    try:
        outfile = h5py.File(filename, 'w')
    except IOError:
        raise IOError('Cannot write file {}'.format(filename))
    
    # Create header...
    try:
        args_githash = ["git", "describe", "--always"]
        githash = subprocess.check_output(args_githash)
    except OSError:
        logging.warning("Git executable not found -- Unable to get Git hash")
        githash = "unknown"
    except subprocess.CalledProcessError: # To be fixed ! Is it because the simulation dir
                                          # may not necessarily be a git repo ?
        logging.warning("Specified command '{}' returned non-zero exit status 128".format(
            " ".join(args_githash)))
        githash = "unknown"
        
        
    outfile.attrs['description'] = \
            'Generated by {} at {} -- Git: {}'.format(getpass.getuser(),
                                                      time.strftime('%Y-%m-%d %H:%M:%S'),
                                                      githash)
    outfile.attrs['metadata'] = metadata
    outfile.attrs['params'] = params
    
    # ... then store attributes and cluster.
    group = outfile.create_group('data')
    for n, c in enumerate(clusters):
        dataset = group.create_dataset('cluster #{}'.format(n), data=c.to_numpyarray(), \
                                       compression='gzip')
        dataset.attrs['description'] = c.metadata
        
    # Close writing session of .hdf file.
    outfile.close()
    
    logging.info('Wrote {} clusters in {}'.format(len(clusters), filename))
