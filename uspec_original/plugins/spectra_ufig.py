#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@student.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.


# System imports
from __future__ import print_function, division, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pf
from wpca import PCA, WPCA, EMPCA

import principal_component_analysis as pri
import uspec.plugins.analyse_uspec as au


class Plugin(object):
    
    def __init__(self, ctx):
        '''
        Constructor
        '''
        self.ctx = ctx

    def reduce(self, ctxList):
        '''
        This function collects the results
        '''

        par = self.ctx.params


        ############
        #  define  #
        ############


        # store spectra and variables in arrays and define other variables
        num = len(ctxList[0].spec_observed)
        quarters_len = int(len(par.survey_lam)-len(par.survey_lam)/4)
        spec_pca = np.zeros((len(ctxList),num))
        spec_nonoise = np.zeros((len(ctxList),num))
        ivar_pca = np.zeros((len(ctxList),num))
        noise = np.zeros((len(ctxList),num))
        level = np.zeros((len(ctxList),num))
        sky = np.zeros((len(ctxList),num))
        sky_true = np.zeros((len(ctxList), num))
        read_noise = np.zeros((len(ctxList), num))
        magnitude_i = np.zeros(len(ctxList))
        z_gals = np.zeros(len(ctxList))

        for idx in range(len(ctxList)):
            spec_pca[idx] = ctxList[idx].spec_observed
            ivar_pca[idx] = ctxList[idx].ivar
            noise[idx] = ctxList[idx].noise_obs
            spec_nonoise[idx] = ctxList[idx].spec_observed_nonoise
            level[idx] = ctxList[idx].spec_level
            sky[idx] = ctxList[idx].sky_erg_obs
            sky_true[idx] = ctxList[idx].sky_obs
            read_noise[idx] = ctxList[idx].readnoise
            magnitude_i[idx] = int(self.ctx.gal_mag_i[idx]*1000.)/1000.
            z_gals[idx] = int(self.ctx.gal_z[idx]*1000.)/1000.

        # save in context for control loop to use
        par.spectra = spec_pca[:,:quarters_len]
        par.lam = par.survey_lam[:quarters_len]

        # select similar USPEC and SDSS galaxies (crude cut in z, mag_i)
        nine_gals = np.zeros((3,len(spec_pca[0][:quarters_len])))
        nine_gals_true = np.zeros((3,len(spec_pca[0][:quarters_len])))
        nine_gals_ivar = np.zeros((3,len(ivar_pca[0][:quarters_len])))
        level_pca = np.zeros((3,len(level[0][:quarters_len])))
        mag_i = np.zeros(3)
        z_nine_gals = np.zeros(3)
        idx_gal = np.zeros(3)
        z_min = [0.475, 0.525, 0.625]
        z_max = [0.5, 0.55, 0.65]
        mag_min = [18.3, 18.5, 19.5]
        mag_max = [18.4, 19.6, 19.6]
        assert_idx = np.zeros(3)

        for idx in range(len(ctxList)):
            for i in range(3):
                if (assert_idx[i] == 0.):
                    if (self.ctx.gal_z[idx] > z_min[i] and self.ctx.gal_z[idx] < z_max[i] and self.ctx.gal_mag_i[idx] < mag_max[i] and self.ctx.gal_mag_i[idx] > mag_min[i]):
                        assert_idx[i] += 1.
                        nine_gals[i] = spec_pca[idx][:quarters_len]
                        nine_gals_true[i] = spec_nonoise[idx][:quarters_len]
                        nine_gals_ivar[i] = ivar_pca[idx][:quarters_len]
                        level_pca[i] = level[idx][:quarters_len]
                        mag_i[i] = int(self.ctx.gal_mag_i[idx]*1000.)/1000.
                        z_nine_gals[i] = int(self.ctx.gal_z[idx]*1000.)/1000.
                        idx_gal[i] = idx


        ############
        #   save   #
        ############


        # outlier rejection
        # if required, put in as weights (WPCA)
        #outliers = np.ones_like(spec_pca)
        #outliers[:,-2:] = 0.
        #outliers[:,:10] = 0.
        #outliers[:,4623] = 0.
        #outliers[:,4173] = 0.


        # conduct PCA and save the results in fits file
        '''
        sample_reconstruction, sample_vecs, sample_coeffs  = pri.pca_sample(WPCA, spec_pca, ncomp=5)
        vecs_pca = np.zeros((len(sample_vecs[0]),len(spec_pca[0][:quarters_len])))
        recon_pca = np.zeros((len(nine_gals),len(spec_pca[0][:quarters_len])))
        recon_all = np.zeros((len(spec_pca),len(spec_pca[0][:quarters_len])))
        for i in range(len(nine_gals)):
            recon_pca[i] = sample_reconstruction[:quarters_len,int(idx_gal[i])]
        for i in range(len(spec_pca)):
            recon_all[i] = sample_reconstruction[:quarters_len, i]
        for idx in range(5):
            vecs_pca[idx] = sample_vecs[:quarters_len,idx]
        hdu_primary = pf.PrimaryHDU()
        cols = pf.ColDefs([
            pf.Column(name=b'spec_true', format=str(len(nine_gals[0])) + 'E', array=spec_nonoise[:,:quarters_len]),
            pf.Column(name=b'spec_true_three', format=str(len(nine_gals[0])) + 'E', array=nine_gals_true),
            pf.Column(name=b'spec', format=str(len(nine_gals[0])) + 'E', array=spec_pca[:,:quarters_len]),
            pf.Column(name=b'spec_three', format=str(len(nine_gals[0])) + 'E', array=nine_gals),
            pf.Column(name=b'ivar', format=str(len(nine_gals[0])) + 'E', array=ivar_pca[:,:quarters_len]),
            pf.Column(name=b'ivar_three', format=str(len(nine_gals[0])) + 'E', array=nine_gals_ivar),
            pf.Column(name=b'spec_pca', format=str(len(nine_gals[0])) + 'E', array=recon_all),
            pf.Column(name=b'spec_pca_three', format=str(len(nine_gals[0])) + 'E', array=recon_pca),
            pf.Column(name=b'p_vecs', format=str(len(nine_gals[0])) + 'E', array=vecs_pca),
            pf.Column(name=b'p_coeffs', format=str(len(sample_coeffs[0])) + 'E', array=sample_coeffs),
            pf.Column(name=b'z', format='E', array=z_gals),
            pf.Column(name=b'z_three', format='E', array=z_nine_gals),
            pf.Column(name=b'mag_i', format='E', array=magnitude_i),
            pf.Column(name=b'mag_i_three', format='E', array=mag_i)])
        hdu_data = pf.BinTableHDU.from_columns(cols)
        del cols
        col1 = pf.Column(name=b'lam', format='E', array=par.survey_lam[:quarters_len])
        cols = pf.ColDefs([col1])
        hdu_lam = pf.BinTableHDU.from_columns(cols)
        hdulist = pf.HDUList([hdu_primary, hdu_data, hdu_lam])
        hdulist.writeto("uspec_compare/uspec.fits", overwrite=True)
        '''

        ###########
        #  PLOTS  #
        ###########


        # conduct PCA
        #pri.plot_results(WPCA, par.survey_lam, spec_pca, len(ctxList), 'uspec', ncomp=5)
        

        # plot spectra with and without noise and the reconstruction including the energy level of the spectra
        #spec_new = {}
        #for idx in range(3):
        #    spec_new['{0}'.format(idx)] = nine_gals_true[idx]
        #    spec_new['{0}'.format(idx+3)] = nine_gals[idx]
        #    spec_new['{0}'.format(idx+6)] = recon_pca[idx]
        #    spec_new['{0}'.format(idx+9)] = level_pca[idx]
        #au.plot_spec(par.survey_lam[:quarters_len], spec_new, vecs_pca, mag_i, z_nine_gals)


        # plot inverse variances
        #au.plot_ivar(par.survey_lam[:quarters_len], nine_gals_ivar, mag_i, z_nine_gals)


        # compare the inverse variance of a random galaxy with its spectrum and plot the inverse variance
        #au.ivar_val(par.survey_lam, spec_pca, spec_nonoise, noise, ivar_pca, len(ctxList))


        # plot kcorrect templates
        #au.plot_kcorrect(par.survey_lam[:4636], par.lam, par.templates_v0)


        # plot nine random uspec cmass spectra
        #au.plot_random_specs(par.survey_lam, spec_pca, spec_nonoise, noise, z_gals, magnitude_i, len(ctxList))


        # plot the read noise and the observed sky spectrum (wo. poisson realization)
        #au.plot_readnoise_sky(par.survey_lam, read_noise, sky, len(ctxList))



    def __str__(self):
        return "reduce and plot plugin"
