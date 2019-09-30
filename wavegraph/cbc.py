# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

# [*]: http://software.ligo.org/docs/lalsuite/lalsimulation/group___l_a_l_sim_inspiral__h.html#gab955e4603c588fe19b39e47870a7b69c

import logging

import lalsimulation
from lalsimulation.lalsimulation import SimInspiralCreateWaveformFlags
from lalsimulation.lalsimulation import GetApproximantFromString
from lalsimulation.lalsimulation import SimInspiralTD

from lal.lal import MSUN_SI as LAL_MSUN_SI   # kg -- mass of the Sun
from lal.lal import PC_SI as LAL_PC_SI       # m -- parsec
from lal.lal import C_SI as LAL_C_SI         # m s^-1 -- speed of light
from lal.lal import G_SI as LAL_G_SI         # m^3 kg^-1 s^-2 -- gravitational constant

LIMIT_ECCENTRICITY = 0.62

PARSER_STR = 'm1={mass1}, m2={mass2}, spin1z={spin1z}, spin2z={spin2z}, ecc={eccentricity}'

class CompactBinaryCoalescence(object):
    """
    A CompactBinaryCoalescence object characterises the gravitational
    wave (GW) chirp signal associated to the coalescence of two
    inspiralling compact objects.
    """
    def __init__(self, mass1, mass2, spin1z, spin2z, eccentricity=0):
        """
        mass1        [float] -- mass of first binary component [Msun]
        mass2        [float] -- mass of second binary component [Msun]
        spin1z       [float] -- Z component of first binary component
        spin2z       [float] -- Z component of second binary component
        eccentricity [float] -- system eccentricity
        """
        self.mass1 = float(mass1)
        self.mass2 = float(mass2)
        self.spin1z = float(spin1z)
        self.spin2z = float(spin2z)
        self.eccentricity = float(eccentricity)
        
    def __attr__(self):
        """ Returns all attributes """
        return self.mass1, self.mass2, self.spin1z, self.spin2z, self.eccentricity
    
    def __str__(self):
        """ Print all attributes """
        return PARSER_STR.format(mass1=self.mass1, mass2=self.mass2,
                                 spin1z=self.spin1z, spin2z=self.spin2z,
                                 eccentricity=self.eccentricity)

    def waveform(self, approximant, freq_min, samp_freq):
        """
        Returns the waveform associated to the binary coalescence. 
        
        Inputs:
        -------
        approximant [str] -- name of the approximant to use. Full list here: [*]
        freq_min  [float] -- lower frequency [Hz]
        samp_freq [float] -- sampling frequency [Hz]

        Outputs:
        --------
        hp [Numpy Array] -- + polarization of the GW [1]
        hc [Numpy Array] -- x polarization of the GW [1]
        """
        # Check inputs.
        circular_nonspinning_approximants = ['TaylorF2', 'SEOBNRv2']
        circular_spinning_approximants = ['SEOBNRv2_ROM_DoubleSpin']
        
        if approximant == 'EccentricTD' and self.eccentricity > LIMIT_ECCENTRICITY:
            logging.error('EccentricTD waveform non valid for eccentricities greater than {}.'.format(LIMIT_ECCENTRICITY))

        if approximant in circular_nonspinning_approximants and \
            self.spin1z != 0 and self.spin2z != 0:
            logging.error('Spins must be zeros. Specified non-spining approximant.')

        # Initialize source and signal parameters
        # See also LAL manual page : http://software.ligo.org/docs/lalsuite/lalsimulation/group___l_a_l_sim_inspiral__c.html
        s1x = s1y = 0                                 # x,y components  of spin vector of object 1
        s2x = s2y = 0                                 # x,y components  of spin vector of object 2
        s1z = self.spin1z                             # z component of spin vector of object 1
        s2z = self.spin2z                             # z component of spin vector of object 2
        D = 1                                         # distance, Mpc
        iota = 0                                      # inclination angle, rad
        phi_ref = 0                                   # gravitational wave phase at end, rad
        longAscNodes = 0                              # longitude of ascending nodes, degenerate with the polarization angle.
        ecc = self.eccentricity                       # eccentricity at reference epoch.
        meanPerAno = 0                                # mean anomaly of periastron.
        deltaT = 1./samp_freq                         # Inverse of sampling frequency, s
        f_min = freq_min                              # starting GW frequency, Hz
        f_ref = 0.                                    # reference GW frequency, Hz
        LALparams = None                              # LAL dictionary containing accessory parameters

        # Initialize approximant
        approximant_flag = GetApproximantFromString(approximant)
        
        # Generate waveform with LALsimulation.
        hp, hc = SimInspiralTD(self.mass1 * LAL_MSUN_SI,
                               self.mass2 * LAL_MSUN_SI,
                               s1x, s1y, s1z, s2x, s2y, s2z,
                               D*LAL_PC_SI*1.0e6,
                               iota,
                               phi_ref,
                               longAscNodes,
                               ecc,
                               meanPerAno,
                               deltaT,
                               f_min,
                               f_ref,
                               LALparams,
                               approximant_flag)
        
        return hp.data.data, hc.data.data

def read_cbc_template_bank(filename):
    """ 
    Read cbc template bank in the .xml or .txt and return a list of 
    CompactBinaryCoalescence objects.
    
    Input:
    ------
    filename -- input .xml or .txt file [string]

    Output:
    ------
    [List of CompactBinaryCoalescence]
    """
    # List of relevant CompactBinaryCoalescence attributes
    attrs = ['mass1', 'mass2', 'spin1z', 'spin2z']

    if not os.path.isfile(filename):
        raise Exception('File {} not found'.format(filename))
    
    logging.info("Reading {}".format(filename))
    
    if os.path.basename(filename).endswith(('.xml', '.xml.gz')):
        
        # Read .xml file
        xmldoc = ligolw_utils.load_filename(filename,
                                            contenthandler=ligolw.LIGOLWContentHandler)
        
        # Read table.
        table = ligolw_table.get_table(xmldoc, \
                                       lsctables.SnglInspiralTable.tableName)
        
        if hasattr(table[0], 'eccentricity'):
            attrs.append('eccentricity')
        
        return [CompactBinaryCoalescence( \
                    *tuple(elem.__getattribute__(attr) for attr in attrs)) \
                    for elem in table]
        
    elif os.path.basename(filename).endswith(('.txt', '.txt.gz')):
        
        data = numpy.loadtxt(filename)

        if data[0].size == 1:
            data = [data]
            
        if data[0].size < 4 or data[0].size > 5:
            raise Exception("Unsupported format -- " \
                            "File must have 4 or 5 fields".format(filename))
        
        if data[0].size == 4:
            logging.info('Input .txt file does not contain eccentricity field.')
        
        if data[0].size == 5:
            logging.info('Input .txt file contains eccentricity field.')
            attrs.append('eccentricity')
        
        return [CompactBinaryCoalescence(*tuple(elem)) for elem in data]
    
    else:
        raise Exception('Unsupported format for CBC template bank')
