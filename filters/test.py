'''
Unit testing
'''
from __future__ import (print_function, division,
                        unicode_literals, absolute_import)
import unittest

import numpy as np
from astropy import units as u
from astropy.analytic_functions import blackbody_nu

from filters import FilterSet

class TestFilterMethods(unittest.TestCase):
    
    def setUp(self):
        self.sdss = FilterSet('sdss')
        T = 5.5e3 * u.K
        wave = np.logspace(2, 5, 100) * u.angstrom
        beam = u.arcsec**2
        flux_density = (blackbody_nu(wave, T) * beam).to(u.Jy)
        self.target = (wave, flux_density)

    def test_center(self):
        centers = np.array([3551, 4686, 6165, 7481, 8931])
        self.assertTrue(np.all(np.isclose(centers, self.sdss.center.to(u.angstrom).value)))
        
    def test_width(self):
        widths = np.array([192.43781296, 388.2933878, 343.66715998, 383.49817568, 525.02721233])
        self.assertTrue(np.all(np.isclose(widths, self.sdss.width.to(u.angstrom).value)))

    def test_flux(self):
        pass
    
if __name__ == '__main__':
    unittest.main()
