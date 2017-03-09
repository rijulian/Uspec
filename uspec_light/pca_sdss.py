#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file fetches the original SDSS CMASS galaxy spectra
#
# System imports
import numpy as np
import astropy.io.fits as aif

# Local modules
from get_plate_id import *


def sdss_data():

    ####################
    #      define      #
    ####################

    # define survey wave length parameters
    wavelength_min                  = 3650.
    wavelength_max                  = 10400.
    num_pixels                      = 4639.

    # input plates
    # original source of the plates is:
    # https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/

    pl_num = np.zeros((10,2))

    #ten random plates
    pl_num[0][0] = int(4288)
    pl_num[0][1] = int(55501)
    pl_num[1][0] = int(4412)
    pl_num[1][1] = int(55912)
    pl_num[2][0] = int(4540)
    pl_num[2][1] = int(55863)
    pl_num[3][0] = int(5044)
    pl_num[3][1] = int(56186)
    pl_num[4][0] = int(5206)
    pl_num[4][1] = int(56033)
    pl_num[5][0] = int(5710)
    pl_num[5][1] = int(56658)
    pl_num[6][0] = int(5962)
    pl_num[6][1] = int(56265)
    pl_num[7][0] = int(6426)
    pl_num[7][1] = int(56334)
    pl_num[8][0] = int(6434)
    pl_num[8][1] = int(56362)
    pl_num[9][0] = int(7111)
    pl_num[9][1] = int(56741)


    # select CMASS galaxies of the plates
    cmass_mask = get_plate_id()


    # define and assign spectra and variables
    sdss_spec_obs = []
    sdss_spec_filtered = []
    spec_obs = []
    lam = []
    observed_spec = [[] for i in range(len(pl_num))]

    # get observed spectra of the CMASS galaxies
    for index in range(len(pl_num)):
        sdss_spec_obs.append('uspec/data/DR13/spPlate-%s-%s.fits'%(int(pl_num[index][0]),int(pl_num[index][1])))
        sdss_spec_filtered.append('uspec/data/DR13/spZbest-%s-%s.fits'%(int(pl_num[index][0]),int(pl_num[index][1])))
        hdulist = aif.open(sdss_spec_obs[index])

        # plate lam
        c0 = hdulist[0].header['coeff0']
        c1 = hdulist[0].header['coeff1']
        num = hdulist[0].header['NAXIS1']
        platelam = np.zeros(num)

        for idx in range(len(platelam)):
            platelam[idx] = c0 + idx * c1
        lam_now = 10**platelam
        lam.append(lam_now)

        # observed spectra
        spec_obs_data = hdulist[0].data[cmass_mask[index]]
        spec_obs.append(spec_obs_data)
        spec_obs[index] = np.array(spec_obs[index])


    # compute survey wave length grid
    coeff0 = np.log10(wavelength_min)
    coeff1 = (np.log10(wavelength_max) - np.log10(wavelength_min)) / num_pixels
    survey_lam = 10 ** (coeff0 + coeff1 * np.arange(num_pixels + 1))

    # interpolate spectra to survey wave length grid
    for idx in range(len(pl_num)):
        for i in range(len(spec_obs[idx])):
            observed_spec[idx].append(np.interp(survey_lam, lam[idx], spec_obs[idx][i]))
        observed_spec[idx] = np.array(observed_spec[idx])

    # merge sdss cmass spectra into one array
    observed_spec_i = observed_spec[0]
    for idx in range(len(pl_num)-1):
        observed_spectrum = np.concatenate((observed_spec_i,observed_spec[idx+1]))
        observed_spec_i = observed_spectrum.tolist()

    # for analysis purposes adjust number of sdss cmass spectra to number of uspec spectra
    uspec_num = 1492
    observed_spectrum = observed_spectrum[:uspec_num]
    quarters_len = len(survey_lam)-len(survey_lam)/4

    # outskirt rejection (get rid of the sky remnant)
    wo_sky_lam = survey_lam[:quarters_len]
    wo_sky_spec = np.zeros((len(observed_spectrum),quarters_len))
    for idx in range(len(observed_spectrum)):
        wo_sky_spec[idx] = observed_spectrum[idx][:quarters_len]


    return survey_lam[:quarters_len], wo_sky_spec
