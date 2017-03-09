#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@student.ethz.ch
#
# according to:
# Copyright (C) 2013 ETH Zurich, lgamper@phys.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.

# System imports
from __future__ import print_function, division, absolute_import, unicode_literals

# External modules
import numpy as np

from ivy.plugin.base_plugin import BasePlugin

class Plugin(BasePlugin):

    def dlam(self, lam):
        """
        Calculate bin widths
        
        :param lam: wavelength grid
        :return dlam: bin widths
        """        

        npix = len(lam)

        dlam = np.zeros_like(lam)
        dlam[:npix-1] = np.abs(lam[1:npix] - lam[:npix-1])
        dlam[npix-1] = np.abs(lam[npix-1] - lam[npix-2])
        

        return dlam

    def ergs17_to_photons(self, lam, exposure_time, telescope_area):
        """
        Compute factor to transform flux in units of 10^-17^ erg/s/cm^2^/Ang to photons per bin
        
        :param lam: wavelength grid in angstrom
        :param exposure_time: exposure time in seconds
        :param telescope_area: telescope area in cm^2

        :return ergs2photonfactor: conversion array, multiply the result with a spectrum in 10^-17^ erg/s/cm^2^/Ang to get number ob photons per bin
        """     

        dlam = self.dlam(lam)

        c = 3.*10**8

        nu = c / (lam*10**-10)

        ergs2photonfactor = 10**(-24) * exposure_time * telescope_area * dlam / (6.616*10**-34 * nu)

        return ergs2photonfactor

    def readnoise(self, lam, instrument_readnoise, instrument_resolution, rn):
        '''
        Read noise for instrument
    
        :param lam: wavelength grid in angstrom
        :param instrument_readnoise: read noise in photon/pixel
        :param instrument_resolution: instrument resolution in angstrom/pixel

        :returns readnoise: read noise in photons
        '''
        dlam = self.dlam(lam)
             
        return instrument_readnoise*dlam*instrument_resolution*np.sqrt(4.*np.pi)*rn
    

    def get_transmission_gal(self, transmission, extinction, fiber_transmission, tr, fiber_switch):
        """
        Transmission for use with the galaxies i.e. with atmospheric extinction.
        
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1 
        :param extinction: atmospheric extinction in 10^-17 erg/s/cm^2^/Ang
        :param fiber_transmission: fiber transmission value
        :param extiction_switch: if False total transmission will be devided by fiber transmission (to compare with SDSS spectra)

        :return transmission: transmission to be used with the galaxies
        """     
        if (fiber_switch == True):
            galtransmission = transmission * np.exp(-extinction/2.5) * tr
        else:
            galtransmission = transmission / fiber_transmission * np.exp(-extinction/2.5) * tr
        
        return galtransmission

    def sky_photons(self, sky, erg_to_photons_factor, transmission, extinction, fiber_transmission, st):
        """
        Creates a Poisson realization of the sky in photons after transmission. 
        
        :param sky: sky spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param erg_to_photons_factor: factor to convert a spectrum on the survey grid from ergs to photons
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1
        :param extinction: atmospheric extinction in 10^-17 erg/s/cm^2^/Ang (if needed)
        :param fiber_transmission: fiber transmission value (if needed)

        :return skyphot: sky in photons after transmission
        """
        skyphot = np.random.poisson(sky * erg_to_photons_factor * transmission * st) #* np.exp(-extinction/2.5)

        return skyphot 

    def instrument_throughput(self, lam, sky, spec, erg_to_photons_factor, transmission, extinction, fiber_transmission, instrument_readnoise,  instrument_resolution, spec_nonoise, rn, tr, st):
        """
        Create a Poisson realization of sky and galaxy photons as they would be measured in the instrument
        
        :param lam: wavelength grid
        :param sky: sky spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param spec: galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang         
        :param erg_to_photons_factor: factor to convert a spectrum on the survey grid from 10^-17^ erg/s/cm^2^/Ang to photons
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1
        :param extinction: atmospheric extinction in 10^-17 erg/s/cm^2^/Ang
        :param fiber_transmission: fiber transmission value
        :param instrument_readnoise: read noise in photons/pixel
        :param instrument_resolution: instrument resolution in angstrom per pixel
        :param spec_nonoise: noiseless galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang

        :return observed: photons of sky and galaxy as they would be measured in the instrument
        :return readnoise: readnoise of the instrument (in photons)
        :return spec_nonoise: noiseless galaxy spectrum as it would be measured in the instrument (still erg)
        """

        # Observed sky emission
        sky_photons = self.sky_photons(sky, erg_to_photons_factor, transmission, extinction, fiber_transmission, st)
        spec_mask = spec > 0
 
        
        # Poisson realisation of galaxy spectrum
        # Spectrum put into poisson realization with a factor of:
        spec_photons = np.zeros_like(spec)
        spec_photons[spec_mask] = np.random.poisson(spec[spec_mask] * erg_to_photons_factor[spec_mask] * self.get_transmission_gal(transmission, extinction, fiber_transmission, tr, fiber_switch=True)[spec_mask])
        
        # Read noise
        readnoise = np.random.normal(size=len(lam)) * self.readnoise(lam, instrument_readnoise,  instrument_resolution, rn)

        observed = spec_photons + sky_photons + readnoise
        
        spec_nonoise = spec_nonoise * self.get_transmission_gal( transmission, extinction, fiber_transmission, tr, fiber_switch=True)

        return observed, readnoise, spec_nonoise

    def throughput_correction(self, sky_photons, photons_after_throughtput, erg_to_photons_factor, transmission, extinction, fiber_transmission, spec_nonoise, tr):
        """
        Correct observed spectrum (in photons) for sky emission and throughput and convert back to 10^-17^ erg/s/cm^2^/Ang

        :param sky_photons: measured sky spectrum in photons
        :param photons_after_throughtput: spectrum (sky + galaxy + read noise) in instrument in photons
        :param erg_to_photons_factor: factor to convert a spectrum on the survey grid from 10^-17^ erg/s/cm^2^/Ang to photons
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1
        :param extinction: atmospheric extinction in 10^-17^ erg/s/cm^2^/Ang
        :param fiber_transmission: fiber transmission value
        :param spec_nonoise: noiseless galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang in instrument

        :return observed: noisy galaxy spectrum after sky substraction and calibration in 10^-17^ erg/s/cm^2^/Ang
        :return spec_nonoise: noiseless galaxy spectrum after calibration in 10^-17^ erg/s/cm^2^/Ang
        """

        factor_mask = erg_to_photons_factor > 0.01
        trans_mask = self.get_transmission_gal(transmission, extinction, fiber_transmission, tr, fiber_switch=False) > 0.0

        observed = photons_after_throughtput - sky_photons
        observed[trans_mask] /= self.get_transmission_gal(transmission, extinction, fiber_transmission, tr, fiber_switch=False)[trans_mask]
        observed[np.logical_not(trans_mask)] = 0.
        observed_ergs = np.zeros_like(observed, dtype="float32")
        observed_ergs[factor_mask] = observed[factor_mask] / erg_to_photons_factor[factor_mask]
        
        spec_nonoise[trans_mask] /= self.get_transmission_gal(transmission, extinction, fiber_transmission, tr, fiber_switch=False)[trans_mask]
        spec_nonoise[np.logical_not(trans_mask)] = 0.

        return observed_ergs, spec_nonoise

    def ivar(self, lam, observed, spec, sky, erg_to_photons_factor, transmission, extinction, fiber_transmission, instrument_readnoise, instrument_resolution, rn, tr, st):
        """
        Generate noise and noise estimations
        
        :param lam: wavelength grid
        :param observed: noisy galaxy spectrum
        :param spec: noise-free spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param sky: sky spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param erg_to_photons_factor: factor to convert a spectrum on the survey grid from 10^-17^ erg/s/cm^2^/Ang to photons
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1
        :param extinction: atmospheric extinction in 10^-17^ erg/s/cm^2^/Ang
        :param fiber_transmission: fiber transmission value                 
        :param instrument_readnoise: read noise in photons/pixel
        :param instrument_resolution: instrument resolution in angstrom per pixel

        :return ivar: inverse variance
        :return noise: noise estimation in 10^-17^ erg/s/cm^2^/Ang
        :return spec_level: energy level of the spectrum in 10^-17^ erg/s/cm^2^/Ang
        """
        #masks
        spec_mask = spec > 0.
        sky_mask = sky > 0.
        
        #photons of spectrum (wo. Poisson realization)
        spec_photons = np.zeros_like(spec)
        spec_photons[spec_mask] = spec[spec_mask] * erg_to_photons_factor[spec_mask] * self.get_transmission_gal(transmission, extinction, fiber_transmission, tr, fiber_switch=True)[spec_mask]

        #sky photons (wo. Poisson realization)
        sky_photons = np.zeros_like(sky)
        sky_photons[sky_mask] = sky[sky_mask] * erg_to_photons_factor[sky_mask] * transmission[sky_mask] * st
        
        #read noise variance
        readnoise_var = np.ones(len(lam)) * abs(self.readnoise(lam, instrument_readnoise,  instrument_resolution, rn))**2
        
        #noise
        n_squared = spec_photons + sky_photons + readnoise_var
        noise = n_squared / erg_to_photons_factor
        
        #fractional inverse variance = signal-to-noise^2
        frac_ivar = spec_photons**2 / n_squared

        #inverse variance
        #polyfit or running mean for ivar level correction
        fit = np.polyfit(lam[spec_mask], spec[spec_mask] * fiber_transmission , 10)
        spec_level = np.polyval(fit, lam[spec_mask])
        ivar = np.zeros_like(observed)
        ivar[spec_mask] = frac_ivar[spec_mask] / spec_level**2 / erg_to_photons_factor**2   #!!!!!!
        
        return ivar, noise, spec_level
        

    def observe(self, lam, spec, sky, extinction, transmission, fiber_transmission, exposure_time, telescope_area, instrument_readnoise, instrument_resolution, spec_nonoise, rn, tr, st):
        """
        Generate observed spectrum and a noise estimation
        
        :param lam: wavelength grid
        :param spec: noise-free spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param sky: sky spectrum in 10^-17^ erg/s/cm^2^/Ang
        :param extinction: atmospheric extinction in 10^-17^ erg/s/cm^2^/Ang
        :param transmission: transmission of atmosphere and instrument in values between 0 and 1
        :param fiber_transmission: fiber transmission value
        :param exposure_time: exposure time in seconds
        :param telescope_area: telescope area in cm^2
        :param instrument_readnoise: read noise in photons/pixel
        :param instrument_resolution: instrument resolution in angstrom/pixel
        :param spec_nonoise: noiseless galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang in instrument

        :return observed: noisy galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang
        :return ivar: inverse variance
        :return sky_erg: sky estimation in 10^-17^ erg/s/cm^2^/Ang
        :return noise: noise estimation in 10^-17^ erg/s/cm^2^/Ang
        :return spec_level: energy level of the spectrum in 10^-17^ erg/s/cm^2^/Ang
        :return readnoise: readnoise of the instrument in photons
        :return spec_nonoise: noiseless galaxy spectrum in 10^-17^ erg/s/cm^2^/Ang
        """

        erg_to_photons_factor = self.ergs17_to_photons(lam, exposure_time, telescope_area)
        factor_mask = erg_to_photons_factor > 0.01

        # Observed photons from sky and galaxy
        photons_after_throughtput, readnoise, spec_nonoise = self.instrument_throughput(lam, sky, spec, erg_to_photons_factor, transmission, extinction, fiber_transmission, instrument_readnoise, instrument_resolution, spec_nonoise, rn, tr, st)

        # Poisson realisation of sky emission used to correct the observed spectrum
        sky_photons = self.sky_photons(sky, erg_to_photons_factor, transmission, extinction, fiber_transmission, st)
            
        # Noise-removed observed spectrum
        observed, spec_nonoise = self.throughput_correction(sky_photons, photons_after_throughtput, erg_to_photons_factor, transmission, extinction, fiber_transmission, spec_nonoise, tr)
        
        # Inverse-variance
        ivar, noise, spec_level = self.ivar(lam, observed, spec, sky, erg_to_photons_factor, transmission, extinction, fiber_transmission, instrument_readnoise, instrument_resolution, rn, tr, st)
        
        sky_erg = sky_photons[factor_mask]/erg_to_photons_factor[factor_mask]
        
        return observed, ivar, sky_erg, noise, spec_level, readnoise, spec_nonoise


    def __call__(self):
        """
        Convert spectra to survey wavelength grid and generate noisy observed spectra
        
        :return spec_observed: noisy galaxy spectrum in survey wavelength grid in units of 10^-17^ erg/s/cm^2^/Ang, corrected for: extinction, throughput, sky substracted        
        :return noise_obs: observed noise in instrument in units of 10^-17^ erg/s/cm^2^/Ang, meaning noise because of: observed sky, read noise, galaxy spectrum
        :return readnoise: read noise of instrument in photons
        :return sky_obs: observed sky eigenspectrum (wo. poisson realization) in units of 10^-17^ erg/s/cm^2^/Ang
        :return sky_erg_obs: observed sky spectrum (incl. poisson realization) in units of 10^-17^ erg/s/cm^2^/Ang
        :return spec_level: energy level of the galaxy spectrum in units of 10^-17^ erg/s/cm^2^/Ang
        :return ivar: inverse variance of the galaxy spectrum
        :return spec_observed_nonoise: noisy galaxy spectrum in survey wavelength grid in units of 10^-17^ erg/s/cm^2^/Ang
        """

        par = self.ctx.params
        
        transmission = np.interp(par.survey_lam, par.lam, self.ctx.params.throughput)

        sky = np.interp(par.survey_lam , par.lam, par.sdss_sky)
        
        extinction = self.ctx.params.survey_airmass * np.interp(par.survey_lam, par.lam, self.ctx.params.extinction)
            
        spec = np.interp(par.survey_lam, self.ctx.lam, self.ctx.spec)
            
        observed, frac_ivar, sky_erg, noise, spec_level, readnoise, spec_nonoise = self.observe(par.survey_lam, spec, sky, extinction, transmission, par.fiber_transmission, par.survey_exposure_time, par.instrument_telescope_area, par.instrument_readnoise, par.instrument_resolution, self.ctx.spec_observed_nonoise, par.rn, par.tr, par.st)


        self.ctx.spec_observed = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.spec_observed = observed 
        
        self.ctx.noise_obs = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.noise_obs = noise
        
        self.ctx.readnoise = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.readnoise = readnoise
        
        self.ctx.sky_obs = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.sky_obs = sky * transmission
        
        self.ctx.sky_erg_obs = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.sky_erg_obs = sky_erg
        
        self.ctx.spec_level = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.spec_level = spec_level
        
        self.ctx.ivar = np.zeros(len(par.survey_lam), dtype="float32")
        self.ctx.ivar = frac_ivar
        
        self.ctx.spec_observed_nonoise = spec_nonoise

        
        del self.ctx.lam
        del self.ctx.spec
        
        print (self.ctx.galaxy_index + 1, "out of", par.numberOfSpectra, "spectra done")

    def __str__(self):
        return "observed spectra module"
