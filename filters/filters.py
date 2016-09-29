'''
Provides classes for modeling photometric filters
'''
from __future__ import (print_function, division,
                        unicode_literals, absolute_import)

import glob

import numpy as np

from scipy.integrate import simps
from scipy.interpolate import interp1d

from astropy import units as u

# root directory of package
from . import _root

class Filter(object):
    '''
    Represents a photometric filter.

    Attributes:
    -----------
    wave : astropy quantity array with units of length
    resp : array of response function (same length as wave)
    
    Methods:
    --------
    __call__ : compute the filter's response to an input spectrum
    '''
    def __init__(self, wave, resp):
        '''
        Parameters
        ----------
        wave : length quantity or array (in which case angstroms are assumed)
        resp : array (equal size to wavelength)
        '''
        try:
            self.wave = wave.to(u.angstrom)
        except AttributeError as e:
            self.wave = wave * u.angstrom
        self.resp = resp
        self.interp = interp1d(self.wave, self.resp,
                               bounds_error=False, fill_value=0)
        # normalization
        norm = simps(resp, wave)
        self.center = simps(resp / norm * wave, wave) * wave.unit # 1st moment
        self.width = np.sqrt(simps(resp / norm * (wave - self.center)**2, wave)) * wave.unit # 2nd moment
        
    def __call__(self, wave, flux):
        '''
        Calculates the integrated flux within the filter response.

        Parameters
        ----------
        wave : wavelength array, astropy quantity of length, or assumed angstroms
        flux : flux array, arbitrary quantity or dimensionless
        '''
        try:
            wave = wave.to(u.angstrom)
        except AttributeError as e:
            wave = wave * u.angstrom
        return simps(self.interp(wave) * flux, wave)


class FilterSet(object):
    '''
    List of filters
    '''
    def __init__(self, filters):
        '''
        filters is either a list of Filter instances or a string with an
        appropriately named filter set (e.g., sdss)
        '''
        if isinstance(filters, str):
            self.filters = []
            if filters.strip().lower() == 'sdss':
                self.name = 'SDSS'
                folder = 'data/sdss/*.dat'
                kwargs = {'usecols': (0, 1)}
            else:
                raise ValueError("Can't find " + filters + " filter set!")
            for filename in glob.glob(_root + folder):
                wave, resp = np.loadtxt(filename, **kwargs).T
                wave = wave * u.angstrom
                self.filters.append(Filter(wave, resp))
        else:
            self.name = ''
            self.filters = filters
        # sort by wave
        idx = np.argsort(self.center.value)
        self.filters = [self.filters[i] for i in idx]

    def __repr__(self):
        return '<' + self.name + ' FilterSet>'
    
    def __getitem__(self, index):
        return self.filters[index]

    def __iter__(self):
        for f in self.filters:
            yield f
            
    def __call__(self, wave, flux):
        return np.array([f(wave, flux) for f in self])
            
    @property
    def center(self):
        return np.array([f.center.value for f in self]) * self[0].wave.unit

    @property
    def width(self):
        return np.array([f.width.value for f in self]) * self[0].wave.unit
