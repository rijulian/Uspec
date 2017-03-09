#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# according to:
# Copyright (C) 2013 ETH Zurich, lgamper@phys.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package
#
# System imports

from __future__ import print_function, division, absolute_import, unicode_literals

import astropy.io.fits as aif
from astropy.io import ascii
import numpy as np

from uspec.plugins.filter_catalog import *

from ivy.plugin.base_plugin import BasePlugin
#from darkskysync.DarkSkySync import DarkSkySync
from clerk.catalog import FitsCatalog
#from clerk.filters import FieldFilter

class Plugin(BasePlugin):
    
    def read_sdss_throughput(self, lam, sdss_throughput_filename):
        '''
        Reads SDSS throughput from fits file and stores the throughput in context.
        Original source of the throughput curves are:
        https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/4412/spFluxcalib-b1-00137140.fits
        https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/4412/spFluxcalib-r1-00137140.fits
        https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/4412/spFrame-b1-00137140.fits
        https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/4412/spFrame-r1-00137140.fits
        For computation see uspec/data/throughput/get_throughput.py
            
        :param lam: wavelength grid in angstrom
        :param sdss_throughput_filename: SDSS spectroscopic filter curves in fits file
            
        :returns throughput: SDSS throughput on lam stored in context
        '''
        hdulist = aif.open(sdss_throughput_filename)
        throughput = hdulist[1].data['throughput'] * 1.66
        loglam = hdulist[1].data['loglam']
        self.ctx.params.survey_airmass = hdulist[0].header['airmass']
        self.ctx.params.throughput = np.interp(lam, 10**(loglam), throughput)
    
    
    def read_sdss_sky(self, lam, sdss_sky_filename):
        '''
        Reads SDSS sky model from fits file and stores the sky model in context.
        Original source of the sky model is:
        https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/
            
        :param lam: wavelength grid in angstrom
        :param sdss_sky_filename: fits file with sky model
            
        :returns sdss_sky: SDSS sky model on lam stored in context
        '''
        hdulist = aif.open(sdss_sky_filename)
        sky_data = hdulist[6].data
        sky = sky_data[0]
        c0 = hdulist[0].header['coeff0']
        c1 = hdulist[0].header['coeff1']
        num = hdulist[6].header['NAXIS1']
        platelam = np.zeros(num)
        for idx in range(len(platelam)):
            platelam[idx] = c0 + idx * c1
        self.ctx.params.sdss_sky = np.interp(lam, 10**platelam, sky)
    
    
    def read_desi_extinction(self, lam, desi_extinction_filename):
        '''
        Reads DESI sky extinction model from ascii file and stores the sky extinction model in context.
        Original source of the sky extinction model:
        https://desi.lbl.gov/trac/browser/code/desimodel/trunk/data/spectra/ZenithExtinction-KPNO.dat

        :param lam: wavelength grid in angstrom
        :param desi_extinction_filename: ascii file with sky extinction model
            
        :returns extinction: DESI sky extinction model on lam stored in context
        '''
        data = ascii.read(desi_extinction_filename)
        extinction = data['EXTINCTION']
        extinctionlam = data['WAVELENGTH']
        self.ctx.params.extinction = np.interp(lam, extinctionlam, extinction)
    
    
    def survey_grid(self, wavelength_min, wavelength_max, num_pixels):
        '''
        Calculates the parameters for a log grid
            
        :param wavelength_min: minimum wavelength fo the grid
        :param wavelength_max: maximum wavelength fo the grid
        :param num_pixels: number of pixels
            
        :returns coeff0: lower end of the survey grid in log stored in context
        :returns coeff1: step size of the survey grid in log stored in context
        :returns npix: number of pixels in survey grid stored in context
        '''
        #create survey grid
        self.ctx.coeff0 = np.log10(wavelength_min)
        self.ctx.coeff1 = (np.log10(wavelength_max) - np.log10(wavelength_min)) / num_pixels
        self.ctx.npix = num_pixels
    
        
    def __call__(self):
        '''
        This module reads in the ufig catalog, the kcorrect templates, the other catalogs (sky etc.) and parameters
        (wave length grid etc.) and creates a wave length grid
        '''
        
        par = self.ctx.params
        
        #dssync = DarkSkySync()
        
        # ufig catalog
        ufig_cat = aif.getdata('uspec/data/catalog/ufig_cmass_cat.fits')
        
        # kcorrect
        kcorrect_template = aif.open('uspec/data/k_nmf_derived.newdefault.fits', memmap=True)

        # get ufig parameters
        mag_r = ufig_cat['mag_r']
        mag_i = ufig_cat['mag_i']
        mag_g = ufig_cat['mag_g']
        coeffs = ufig_cat['coeffs']
        z = ufig_cat['z']
        
        # cmass cut
        mask = CmassFilter(mag_g, mag_r, mag_i, mag_i, mag_i, z)
        mask_ufig = mask.get_cmass()
        
        # save data to context
        self.ctx.gal_mag_r = mag_r[mask_ufig]
        self.ctx.gal_mag_i = mag_i[mask_ufig]
        self.ctx.gal_coeffs = coeffs[mask_ufig]
        self.ctx.gal_z = z[mask_ufig]
        
        # wavelength grid
        par.lam = kcorrect_template[11].data
        
        # templates
        par.templates_v0 = np.array(kcorrect_template[5].data, dtype='float32')
        
        # spectroscopic transmission curve sdss
        par.sdss_throughput_filename = 'uspec/data/throughput/eff_sdss.fits'
        self.read_sdss_throughput(par.lam, par.sdss_throughput_filename)
        
        # sdss sky
        sdss_sky_filename = 'uspec/data/DR13/spPlate-4412-55912.fits'
        self.read_sdss_sky(par.lam, sdss_sky_filename)
        
        # extinction
        desi_extinction_filename = 'uspec/data/ZenithExtinction-KPNO.dat'
        self.read_desi_extinction(par.lam, desi_extinction_filename)

        # survey parameters
        par.instrument_wavelength_range = [3650., 10400.]    #boss
        #par.instrument_wavelength_range = [3800., 9200.]    #sdss
        
        par.survey_exposure_time            = 2700.     #sec per plate (approx 45 min)
        par.instrument_num_pixels 	        = 4649.     #see sdss files
        par.instrument_readnoise		    = 3.        #photons / pixel
        par.instrument_resolution           = (par.instrument_wavelength_range[1]-par.instrument_wavelength_range[0])/(par.instrument_num_pixels+1.)  #  angstrom / pixel 
        par.instrument_telescope_area       = np.pi * (125.**2 - 65.0**2)   # cm^2
        par.fiber_transmission              = 0.537     #0.47 (modelmag) #0.434 (dered) #0.537 (petro)  #see pca_sdss.py
                                                        #comparison with sdss: 0.48903208971

        par.sky_model = 'boss'
        #par.sky_model = 'desi'


        #control loop parameters
        par.rn = 1.2391  #abs(par.theta[0])          #readnoise control constant multiply in readnoise
        par.tr = 0.9364  #abs(par.theta[1])          #throughput control constant (maybe shape) multiply in get_transmission_gal
        par.st = 1.2152  #abs(par.theta[2])          #throughput control constant (maybe shape) multiply in sky_photons


        self.survey_grid(par.instrument_wavelength_range[0], par.instrument_wavelength_range[1], par.instrument_num_pixels)
    
        survey_lam = 10 ** (self.ctx.coeff0 + self.ctx.coeff1 * np.arange(self.ctx.npix + 1))
        survey_mask = np.interp(survey_lam, par.lam, self.ctx.params.throughput) > 5e-4
        par.survey_lam = survey_lam[survey_mask]


    def __str__(self):
        return "import ufig data"