'''
Provides classes for modeling photometric filters
'''
from __future__ import (print_function, division,
                        unicode_literals, absolute_import)

import os
import glob

import numpy as np

from scipy.integrate import simps
from scipy.interpolate import interp1d

from astropy import units as u

# root directory of package
_root = os.path.abspath(os.path.dirname(__file__)) + '/'

# AB flux zeropoint
zp_AB = 3631 * u.Jy

class Filter(object):
    '''
    Represents a photometric filter.

    Attributes:
    -----------
    wave : quantity array with units of length
    resp : array of response function (same length as wave)
    name : str, assigned name
    zp : quantity, flux zeropoint 

    Methods:
    --------
    __call__ : compute the filter's response to an input spectrum
    interp : regrid onto a new spectral domain
    '''
    def __init__(self, wave, resp, name='', zp=zp_AB):
        '''
        Parameters
        ----------
        wave : length quantity or array (in which case angstroms are assumed)
        resp : array (equal size to wavelength) of total quantum efficiency
        '''
        try:
            self.wave = wave.to(u.angstrom)
        except AttributeError as e:
            wave = wave * u.angstrom
        self.wave = wave
        self.resp = resp
        assert np.all(resp <= 1)
        self.name = name
        self.zp = zp
        self.norm = simps(resp, wave)
        # 1st moment        
        self.center = simps(resp / self.norm * wave, wave) * wave.unit
        # 2nd moment        
        self.width = np.sqrt(simps(resp / self.norm * (wave - self.center)**2, wave)) * wave.unit 

    def __repr__(self):
        return '<Filter: ' + self.name + '>'
    
    def __call__(self, wave, flux_density, mag=False):
        '''
        Calculates the integrated flux within the filter response.

        Parameters
        ----------
        wave : wavelength array, astropy quantity of length, or assumed angstroms
        flux_density : flux density array, quantity of flux density, or assumed janskys
        mag : bool, if true, use the zeropoint to assign a magnitude instead of returning flux

        Returns
        -------
        flux : quantity of flux (or magnitude if selected)
        '''
        try:
            freq = wave.to(u.Hz, equivalencies=u.spectral())
        except AttributeError as e:
            wave = wave * u.angstrom
            freq = wave.to(u.Hz, equivalencies=u.spectral())
        try:
            flux_density = flux_density.to(u.Jy, equivalencies=u.spectral_density(wave))
        except AttributeError as e:
            flux_density = flux_density * u.Jy
        # sort by increasing frequency
        idx = np.argsort(freq)
        new_resp = self.interp(wave)[idx]
        freq = freq[idx]
        flux_density = flux_density[idx]
        flux = simps(new_resp * flux_density, freq) * flux_density.unit * freq.unit
        if mag:
            return -2.5 * np.log10((flux.to(self.zp.unit * u.Hz) /
                                    (simps(self.zp * new_resp, freq) *
                                     self.zp.unit * freq.unit)).to(u.dimensionless_unscaled))
        return flux
                
    def interp(self, new_wave):
        '''
        Interpolate the response curve onto a new spectral grid.

        Parameters
        ----------
        new_wave : length quantity array, or assumed angstroms

        Returns
        -------
        new_resp : array
        '''
        f = interp1d(self.wave, self.resp, bounds_error=False, fill_value=0)
        try:
            new_wave = new_wave.to(u.angstrom, equivalencies=u.spectral())
        except AttributeError as e:
            new_wave = new_wave * u.angstrom
        return f(new_wave)
            
        
class FilterSet(object):
    '''
    List of filters
    '''
    def __init__(self, filters, name=''):
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
            elif filters.strip().lower() in ['cfht', 'megacam', 'cfht/megacam']:
                folder = 'data/cfht/megacam/*.dat'
                kwargs = {'usecols': (0, 1), 'delimiter': ','}
                self.name = 'CFHT/MegaCam'
            else:
                raise ValueError("Can't find " + filters + " filter set!")
            for filename in glob.glob(_root + folder):
                print(filename)
                wave, resp = np.loadtxt(filename, **kwargs).T
                wave = wave * u.angstrom
                self.filters.append(Filter(wave, resp))
        else:
            self.name = name
            self.filters = filters
        assert len(filters) > 0
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
            
    def __call__(self, wave, flux, mag=False):
        return np.array([f(wave, flux, mag) for f in self])
            
    @property
    def center(self):
        return np.array([f.center.value for f in self]) * self[0].wave.unit

    @property
    def width(self):
        return np.array([f.width.value for f in self]) * self[0].wave.unit

