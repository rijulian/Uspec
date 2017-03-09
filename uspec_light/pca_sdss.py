#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file conducts a PCA of the SDSS spectra stored in uspec/data/DR13
#
# System imports
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as aif
from wpca import PCA, WPCA, EMPCA

# Local modules
from get_plate_id import *
import principal_component_analysis as pri
import analyse_sdss as asd


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

    '''
    #add another eight plates, if needed
    pl_num[10][0] = int(3760)
    pl_num[10][1] = int(55268)
    pl_num[11][0] = int(4034)
    pl_num[11][1] = int(55635)
    pl_num[12][0] = int(4998)
    pl_num[12][1] = int(55722)
    pl_num[13][0] = int(5355)
    pl_num[13][1] = int(56009)
    pl_num[14][0] = int(5899)
    pl_num[14][1] = int(56038)
    pl_num[15][0] = int(6303)
    pl_num[15][1] = int(56539)
    pl_num[16][0] = int(6618)
    pl_num[16][1] = int(56401)
    pl_num[17][0] = int(6837)
    pl_num[17][1] = int(56442)
    '''

    '''
    # if need be redshift the spectra back to local system
    def redshift_back(lam, spec, z):
        """
        Redshift spectrum back to z=0 (local)

        :param lam: wavelength grid, log or linear
        :param spec: spectrum
        :param z: redshift
        :returns lam_shifted: back redshifted wavelength grid
        :returns spectrum_shifted: galaxy spectrum back redshifted
        """
        # take spectrum and redshift it
        # spectrum has units originally (erg/s/cm^2/A)
        d_lam = lam[1] - lam[0]            # size of bin in A
        lam_min = lam - 0.5 * d_lam     # edges of bin in A
        lam_max = lam + 0.5 * d_lam     # edges of bin in A
        lam_min_shift = lam_min / (1 + z)                    # redshifted edges of bin in A
        lam_max_shift = lam_max / (1 + z)                    # redshifted edges of bin in A
        lam_shift = (lam_min_shift + lam_max_shift) / 2.0     # redshifted mean of bin in A
        spec_shift = spec * ( 1 + z )
        return lam_shift,  spec_shift
    '''


    # select CMASS galaxies of the plates
    cmass_mask, z_plate, mag_plate, magr_plate = get_plate_id()


    # define and assign spectra and variables
    sdss_spec_obs = []
    sdss_spec_filtered = []
    spec_obs = []
    spec_eigen = []
    ivar = []
    lam = []
    observed_spec = [[] for i in range(len(pl_num))]
    observed_eigen = [[] for i in range(len(pl_num))]
    observed_ivar = [[] for i in range(len(pl_num))]

    # get observed spectra of the CMASS galaxies
    for index in range(len(pl_num)):
        sdss_spec_obs.append('uspec/data/DR13/spPlate-%s-%s.fits'%(int(pl_num[index][0]),int(pl_num[index][1])))
        sdss_spec_filtered.append('uspec/data/DR13/spZbest-%s-%s.fits'%(int(pl_num[index][0]),int(pl_num[index][1])))
        hdulist = aif.open(sdss_spec_obs[index])
        hdulist2 = aif.open(sdss_spec_filtered[index])

        # plate lam
        c0 = hdulist[0].header['coeff0']
        c1 = hdulist[0].header['coeff1']
        num = hdulist[0].header['NAXIS1']
        platelam = np.zeros(num)

        for idx in range(len(platelam)):
            platelam[idx] = c0 + idx * c1
        lam_now = 10**platelam

        # observed spectra
        spec_obs_data = hdulist[0].data[cmass_mask[index]]
        spec_obs.append(spec_obs_data)
        spec_obs[index] = np.array(spec_obs[index])

        # eigen spectra
        spec_eigen_data = hdulist2[2].data[cmass_mask[index]]
        spec_eigen.append(spec_eigen_data)
        spec_eigen[index] = np.array(spec_eigen[index])

        # inverse variance of the spectra
        ivar_data = hdulist[1].data[cmass_mask[index]]
        ivar.append(ivar_data)
        ivar[index] = np.array(ivar[index])

        ## back-redshift the spectra, if needed
        #for i in range(len(spec_obs[index])):
        #    lam_inter, spec_obs[index][i] = redshift_back(lam_now, spec_obs[index][i], z_plate[index][i])
        #    lam[index].append(lam_inter)
        lam.append(lam_now)


    # compute survey wave length grid
    coeff0 = np.log10(wavelength_min)
    coeff1 = (np.log10(wavelength_max) - np.log10(wavelength_min)) / num_pixels
    survey_lam = 10 ** (coeff0 + coeff1 * np.arange(num_pixels + 1))

    # interpolate spectra to survey wave length grid
    for idx in range(len(pl_num)):
        for i in range(len(spec_obs[idx])):
            observed_spec[idx].append(np.interp(survey_lam, lam[idx], spec_obs[idx][i]))
            observed_eigen[idx].append(np.interp(survey_lam, lam[idx], spec_eigen[idx][i]))
            observed_ivar[idx].append(np.interp(survey_lam, lam[idx], ivar[idx][i]))
        observed_spec[idx] = np.array(observed_spec[idx])
        observed_eigen[idx] = np.array(observed_eigen[idx])
        observed_ivar[idx] = np.array(observed_ivar[idx])


    # merge sdss cmass spectra into one array
    observed_spec_i = observed_spec[0]
    observed_eigen_i = observed_eigen[0]
    observed_ivar_i = observed_ivar[0]

    for idx in range(len(pl_num)-1):
        observed_spectrum = np.concatenate((observed_spec_i,observed_spec[idx+1]))
        observed_eigenspec = np.concatenate((observed_eigen_i,observed_eigen[idx+1]))
        observed_total_ivar = np.concatenate((observed_ivar_i,observed_ivar[idx+1]))
        observed_spec_i = observed_spectrum.tolist()
        observed_eigen_i = observed_eigenspec.tolist()
        observed_ivar_i = observed_total_ivar.tolist()


    # for analysis purposes adjust number of sdss cmass spectra to number of uspec spectra
    uspec_num = 1492
    observed_spectrum = observed_spectrum[:uspec_num]
    observed_eigenspec = observed_eigenspec[:uspec_num]
    observed_total_ivar = observed_total_ivar[:uspec_num]

    # merge sdss cmass parameters into one array
    number = len(observed_spectrum)
    z_plate = np.concatenate(z_plate)
    mag_plate = np.concatenate(mag_plate)
    magr_plate = np.concatenate(magr_plate)
    quarters_len = len(survey_lam)-len(survey_lam)/4


    # select similar USPEC and SDSS galaxies (crude cut in z, mag_i)
    nine_gals = np.zeros((3,len(observed_spectrum[0][:quarters_len])))
    nine_gals_true = np.zeros((3,len(observed_eigenspec[0][:quarters_len])))
    nine_gals_ivar = np.zeros((3,len(observed_total_ivar[0][:quarters_len])))
    mag_i = np.zeros(3)
    mag_r = np.zeros(3)
    z_nine_gals = np.zeros(3)
    idx_gal = np.zeros(3)
    z_min = [0.475, 0.525, 0.625]
    z_max = [0.5, 0.55, 0.65]
    mag_min = [18.3, 18.5, 19.5]
    mag_max = [18.5, 19.6, 19.6]
    assert_idx = np.zeros(3)

    for idx in range(number):
        for i in range(3):
            if (assert_idx[i] == 0.):
                if (z_plate[idx] > z_min[i] and z_plate[idx] < z_max[i] and mag_plate[idx] < mag_max[i] and mag_plate[idx] > mag_min[i]):
                    assert_idx[i] = 1.
                    nine_gals[i] = observed_spectrum[idx][:quarters_len]
                    nine_gals_true[i] = observed_eigenspec[idx][:quarters_len]
                    nine_gals_ivar[i] = observed_total_ivar[idx][:quarters_len]
                    mag_i[i] = int(mag_plate[idx]*1000.)/1000.
                    mag_r[i] = int(magr_plate[idx]*1000.)/1000.
                    z_nine_gals[i] = int(z_plate[idx]*1000.)/1000.
                    idx_gal[i] = idx


    # outskirt rejection (get rid of the sky remnant)
    wo_sky_lam = survey_lam[:quarters_len]
    wo_sky_spec = np.zeros((len(observed_spectrum),quarters_len))
    wo_sky_true = np.zeros((len(observed_eigenspec),quarters_len))
    wo_sky_ivar = np.zeros((len(observed_total_ivar),quarters_len))
    for idx in range(len(observed_spectrum)):
        wo_sky_spec[idx] = observed_spectrum[idx][:quarters_len]
        wo_sky_true[idx] = observed_eigenspec[idx][:quarters_len]
        wo_sky_ivar[idx] = observed_total_ivar[idx][:quarters_len]


    '''
    # estimate r-band magnitude as a comparison, if need be
    # calculate average fiber transmission
    fiber_transmission_r = np.zeros(number)
    for idx in range(number):
        modelmag_r = asd.check_mag_r(survey_lam, observed_eigenspec[idx]/1e17, magr_plate[idx])
        fiber_transmission_r[idx] = 10**(2./5.*(magr_plate[idx]-modelmag_r))
    fiber_transmission_r = np.mean(fiber_transmission_r)
    print(fiber_transmission_r)
    '''



    ########################
    #        Plots         #
    ########################


    # conduct PCA
    #pri.plot_results(WPCA, wo_sky_lam, wo_sky_spec, len(wo_sky_spec), 'sdss', ncomp=5)


    # sample galaxy plot
    #asd.sample_plot(wo_sky_lam, wo_sky_spec, len(wo_sky_spec))


    # plot random sdss cmass spectra
    #asd.plot_random_specs(wo_sky_lam, wo_sky_spec, wo_sky_true, z_plate, mag_plate, len(wo_sky_spec))


    # plot inverse variances
    #asd.plot_ivar(wo_sky_lam, nine_gals_ivar, mag_i, z_nine_gals)


    # plot several different pca possibilities
    #Y_wosky_ivar, vecs_wosky_ivar, coeffs_wosky_ivar = pri.pca_sample(WPCA, wo_sky_lam, wo_sky_spec, weights=wo_sky_ivar, ncomp=5)
    #Y_wosky, vecs_wosky, coeffs_wosky = pri.pca_sample(WPCA, wo_sky_lam, wo_sky_spec, ncomp=5)
    #Y_normal, vecs_normal, coeffs_normal = pri.pca_sample(WPCA, survey_lam, observed_spectrum, ncomp=5)
    #Y_ivar, vecs_ivar, coeffs_ivar = pri.pca_sample(WPCA, survey_lam, observed_spectrum, weights=observed_total_ivar, ncomp=5)


    # plot the four different PCA possibilities
    #asd.plot_reconstructed(survey_lam, wo_sky_lam, Y_normal, Y_ivar, Y_wosky, Y_wosky_ivar, number)
    #asd.plot_vectors(survey_lam, wo_sky_lam, vecs_normal, vecs_ivar, vecs_wosky, vecs_wosky_ivar, number)



    ##################
    #      save      #
    ##################

    '''
    # conduct PCA and save in fits file
    sample_reconstruction, sample_vecs, sample_coeffs  = pri.pca_sample(WPCA, wo_sky_spec, ncomp=5)
    pca_rec = np.zeros((3,len(observed_spectrum[0][:quarters_len])))
    pca_rec_all = np.zeros((len(observed_spectrum),len(observed_spectrum[0][:quarters_len])))
    vecs = np.zeros((len(sample_vecs[0]),len(observed_spectrum[0][:quarters_len])))
    for i in range(3):
        pca_rec[i] = sample_reconstruction[:,int(idx_gal[i])]
    for i in range(len(wo_sky_spec)):
        pca_rec_all[i] = sample_reconstruction[:,i]
    for idx in range(5):
        vecs[idx] = sample_vecs[:,idx]

    hdu_primary = aif.PrimaryHDU()
    cols = aif.ColDefs([
        aif.Column(name=b'spec_true', format=str(len(wo_sky_true[0])) + 'E', array=wo_sky_true),
        aif.Column(name=b'spec_true_three', format=str(len(nine_gals_true[0])) + 'E', array=nine_gals_true),
        aif.Column(name=b'spec', format=str(len(wo_sky_spec[0])) + 'E', array=wo_sky_spec),
        aif.Column(name=b'spec_three', format=str(len(nine_gals[0])) + 'E', disp=bytes(b'F10.6'), array=nine_gals),
        aif.Column(name=b'ivar', format=str(len(wo_sky_ivar[0])) + 'E', disp=bytes(b'F10.6'), array=wo_sky_ivar),
        aif.Column(name=b'ivar_three', format=str(len(nine_gals_ivar[0])) + 'E', disp=bytes(b'F10.6'), array=nine_gals_ivar),
        aif.Column(name=b'spec_pca', format=str(len(pca_rec_all[0])) + 'E', disp=bytes(b'F10.6'), array=pca_rec_all),
        aif.Column(name=b'spec_pca_three', format=str(len(pca_rec[0])) + 'E', disp=bytes(b'F10.6'), array=pca_rec),
        aif.Column(name=b'p_vecs', format=str(len(vecs[0])) + 'E', disp=bytes(b'F10.6'), array=vecs),
        aif.Column(name=b'p_coeffs', format=str(len(sample_coeffs[0])) + 'E', disp=bytes(b'F10.6'), array=sample_coeffs),
        aif.Column(name=b'z', format='E', disp=bytes(b'F10.6'), array=z_plate[:number]),
        aif.Column(name=b'z_three', format='E', disp=bytes(b'F10.6'), array=z_nine_gals),
        aif.Column(name=b'mag_i', format='E', disp=bytes(b'F10.6'), array=mag_plate[:number]),
        aif.Column(name=b'mag_i_three', format='E', disp=bytes(b'F10.6'), array=mag_i)])
    hdu_data = aif.BinTableHDU.from_columns(cols)
    del cols
    col1 = aif.Column(name=b'lam', format='E', array=survey_lam[:quarters_len])
    cols = aif.ColDefs([col1])
    hdu_lam = aif.BinTableHDU.from_columns(cols)
    hdulist = aif.HDUList([hdu_primary, hdu_data, hdu_lam])
    hdulist.writeto("uspec_compare/sdss.fits", overwrite=True)
    '''

    return survey_lam[:quarters_len], wo_sky_spec
