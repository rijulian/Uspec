#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# according to:
# Copyright (C) 2013 ETH Zurich, lgamper@phys.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# This file creates a redshifted spectrum


# System imports
from __future__ import print_function, division, absolute_import, unicode_literals

# External modules
import numpy as np

from ivy.plugin.base_plugin import BasePlugin

class Plugin(BasePlugin):
    
    def redshift_spectrum(self,lam, spec, z):
        """
        Redshift spectrum

        :param lam: wavelength grid, log or linear
        :param spec: spectrum
        :param z: redshift

        :returns lam_shifted: redshifted wavelength grid
        :returns spectrum_shifted: galaxy spectrum redshifted to galaxy redshift
        """
        
        # take spectrum and redshift it
        # spectrum has units originally (erg/s/cm^2/A)
        d_lam = lam[1] - lam[0]            # size of bin in A
        lam_min = lam - 0.5 * d_lam        # edges of bin in A
        lam_max = lam + 0.5 * d_lam        # edges of bin in A
        
        lam_min_shift = lam_min * (1.+ z)                     # redshifted edges of bin in A
        lam_max_shift = lam_max * (1.+ z)                     # redshifted edges of bin in A
        lam_shifted = (lam_min_shift + lam_max_shift) / 2.0     # redshifted mean of bin in A
        
        spec_redshifted = spec / (1.+ z)
        
        return lam_shifted,  spec_redshifted
    
    def redshifted_spectrum(self, lam, templates, coeffs, z):
        """
        Creates a redshifted spectrum using templates, coeffs and a redshift

        :param lam: wavelength grid
        :param templates: templates to be used with coeffs
        :param coeffs: coefficients for the templates
        :param z: redshift

        :return lam_redshifted: wavelength grid of redshifted spectrum
        :return spectrum_reshifted: redshifted spectrum
        """

        #take first five coeffs and produce redshifted spectrum
        #spec = np.array(np.matrix(coeffs)*np.matrix(templates))[0]

        spec = np.zeros_like(templates[0])

        for i in xrange(len(coeffs)):
            spec += coeffs[i] * templates[i]

        lam_redshifted, spectrum_reshifted = self.redshift_spectrum(lam, spec, z)
        
        return lam_redshifted, spectrum_reshifted
    
    def calculate_magnitude_r(self,lam, spec):
        """
        Calculate the r-band magnitude of a given spectrum for the DES camera

        :param lam: wavelength grid, log or linear
        :param spec: spectrum

        :return mag_r: r magnitude for the DES camera
        """

        c = 3 * 1e8                                              # m/s
    
        # given filter specification, total throughput and a spectrum, calcualte magnitude
        # nu lam = c, nu = c / lam, d_nu = c / lam^2 d_lam
        #TODO: find the right values for sdss or use the filter curves

        l_mean = 6250.  #sdss r-band mean
        delta_l = 1400. #sdss r-band FWHM
        l_min = l_mean - 0.5 * delta_l
        l_max = l_mean + 0.5 * delta_l
    
        mask = (lam > l_min)*(lam < l_max)
        lam = lam[mask]
        d_lam = lam[1:] - lam[:-1]           # A
        lam = (lam[1:] + lam[:-1]) / 2       # A
        spec = spec[mask]
        spec = (spec[1:] + spec[:-1]) / 2    # erg/s/cm^2/A

        #https://en.wikipedia.org/wiki/AB_magnitude
        int1 = np.sum(spec * lam / (c*1e10) * d_lam )
        int2 = np.sum(1./lam*d_lam)  # 1/s=Hz
        # write this again
        mag_r = -2.5 * np.log10( int1 / int2 ) - 48.6
    
        return mag_r

    def check_mag_r(self, lam, spec, mag_r):
        """
        Calculate the r band magnitude and compares it to given magnitude, prints a warning if they differ by more than 5%

        :param lam: wavelength grid
        :param spec: spectrum
        :param mag_r: magnitude to compare with

        :return mag_r_calc: calculated mag_r of the spectrum
        """

        mag_r_calc = self.calculate_magnitude_r(lam, spec)
        if np.abs(mag_r - mag_r_calc) / mag_r > 0.05:
            print("mag_r a bit off (> 5%): ", np.abs(mag_r - mag_r_calc) / mag_r)
        
        return mag_r_calc
    
    def redshift_kcorrect_spectrum(self):
        """
        Create a redshifted using data from context

        :return lam: wavelength grid of redshifted ufig spectrum stored in context
        :return spec: redshifted ufig spectrum stored in context
        :return mag_r_calc: calculated r-band magnitude stored in context
        """
        
        par = self.ctx.params
        
        # *1e17 to be compatible with sdss
        self.ctx.lam, self.ctx.spec = self.redshifted_spectrum(par.lam, par.templates_v0, self.ctx.coeffs * 1e17, self.ctx.z)
        self.ctx.mag_r_calc = self.check_mag_r(self.ctx.lam, self.ctx.spec / 1e17, self.ctx.mag_r)


    def __call__(self):
        """
        Generate a redshifted spectrum (10^{17} fold too large)
        """
            
        par = self.ctx.params
            
        self.redshift_kcorrect_spectrum()
        
        #store noiseless spectrum in context
        self.ctx.spec_observed_nonoise = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.spec_observed_nonoise = np.interp(par.survey_lam, self.ctx.lam, self.ctx.spec)


    def __str__(self):
        return "simulate spectrum"