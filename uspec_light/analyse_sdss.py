#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file contains analysing tools for SDSS spectra
#
# System imports
import numpy as np
import matplotlib.pyplot as plt


def sample_plot(survey_lam, observed_spectrum, galnum):
    '''
    Plot a random sample galaxy spectrum

    :param survey_lam: wave length grid
    :param observed_spectrum: observed galaxy spectra
    :param galnum: number of galaxies available
    '''
    ran_num = np.random.choice(galnum)
    plt.figure()
    plt.title('Random SDSS galaxy spectrum')
    plt.plot(survey_lam, observed_spectrum[ran_num], color="k", linewidth=0.5)
    plt.xlabel('$\lambda \ [\AA]$')
    plt.ylim([-2, 5])
    plt.savefig("plots/sdss/sample_gal.pdf")
    plt.show()


def plot_random_specs(survey_lam, spec, spec_nonoise, z_gals, mag_i, num_gals):
    '''
    Plot random noisy and noiseless galaxy spectra

    :param survey_lam: wave length grid
    :param spec: galaxy spectra
    :param spec_nonoise: noiseless galaxy spectra
    :param z_gals: redshifts
    :param mag_i: i-band magnitudes
    :param num_gals: number of galaxies available
    '''
    ran_idx = np.random.choice(num_gals, 9)
    plt.figure(figsize=(25,10))
    idxs = 1
    for idx in ran_idx:
        plt.subplot(3,3,idxs)
        plt.title(r'$z$ = %s , $mag_i$ = %s'%(z_gals[idx], mag_i[idx]))
        plt.plot(survey_lam, spec[idx], color="k", zorder=2, linewidth=0.5, label="Noisy spectrum of galaxy")
        plt.plot(survey_lam, spec_nonoise[idx], color="r", zorder=3, linewidth=1., label="Noise-free spectrum of galaxy")
        plt.ylim([-2.,10.])
        plt.legend(loc='best',prop={'size':6})
        plt.ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
        if (idxs > 6):
            plt.xlabel('$\lambda \ [\AA]$')
        idxs += 1
    plt.suptitle('Random SDSS CMASS galaxies', fontsize=18)
    plt.savefig("plots/sdss/spec.pdf")
    plt.show()


def plot_ivar(survey_lam, nine_gals_ivar, mag_i, z_nine_gals):
    '''
    Plot inverse variances of three predefined galaxies

    :param survey_lam: wave length grid
    :param nine_gals_ivar: inverse variances
    :param mag_i: i-band magnitude
    :param z_nine_gals: redshifts
    '''
    fig, a = plt.subplots(1, 3, figsize=(20, 5))
    a = a.ravel()
    for idx,ax in enumerate(a):
        ax.plot(survey_lam, nine_gals_ivar[idx], c='black', lw=0.5)
        ax.text(3600, 25, r'$mag_i$: %s, $z$: %s' % (mag_i[idx], z_nine_gals[idx]))
        ax.set_xlabel('$\lambda \ [\AA]$')
        ax.set_ylim([0, 30])
    fig.suptitle('Inverse Variances', fontsize=16)
    plt.savefig("plots/sdss/multiple_ivar.pdf")
    plt.show()


def plot_reconstructed(x_lam1, x_lam2, y_recon1, y_recon2, y_recon3, y_recon4, number):
    '''
    Plot reconstructed spectrum of a random galaxy done with a PCA in four different ways

    :param x_lam1: long wave length grid
    :param x_lam2: short wave length grid
    :param y_recon1: reconstructed spectrum done with sky, without weights
    :param y_recon2: reconstructed spectrum done with sky, with weights
    :param y_recon3: reconstructed spectrum done without sky, without weights
    :param y_recon4: reconstructed spectrum done without sky, with weights
    :param number: number of available galaxies
    '''
    ran_num = np.random.choice(number)
    fig, ax = plt.subplots(2, 2, figsize=(16, 6))
    ax[0, 0].plot(x_lam1, y_recon1[:, ran_num], c='black', lw=0.5)
    ax[0, 0].set_ylim([-2, 5])
    ax[1, 0].plot(x_lam1, y_recon2[:, ran_num], c='black', lw=0.5)
    ax[1, 0].set_ylim([-2, 5])
    ax[0, 1].plot(x_lam2, y_recon3[:, ran_num], c='black', lw=0.5)
    ax[0, 1].set_ylim([-2, 5])
    ax[1, 1].plot(x_lam2, y_recon4[:, ran_num], c='black', lw=0.5)
    ax[1, 1].set_ylim([-2, 5])

    ax[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0, 1].xaxis.set_major_formatter(plt.NullFormatter())

    ax[0, 0].set_title('PCA with sky, no weights')
    ax[1, 0].set_title('PCA with sky, ivar as weights')
    ax[0, 1].set_title('PCA without sky, no weights')
    ax[1, 1].set_title('PCA without sky, ivar as weights')

    fig.suptitle('Reconstructed spectrum of random galaxy (PCA with %s gals)' % number, fontsize=16)
    plt.savefig("plots/sdss/PCA_SDSS_reconstructed.pdf")
    plt.show()


def plot_vectors(x_lam1, x_lam2, y_vec1, y_vec2, y_vec3, y_vec4, number):
    '''
    Plot first five principal vectors of PCA done in four different ways

    :param x_lam1: long wave length grid
    :param x_lam2: short wave length grid
    :param y_recon1: reconstructed spectrum done with sky, without weights
    :param y_recon2: reconstructed spectrum done with sky, with weights
    :param y_recon3: reconstructed spectrum done without sky, without weights
    :param y_recon4: reconstructed spectrum done without sky, with weights
    :param number: number of available galaxies
    '''
    fig, ax = plt.subplots(2, 2, figsize=(16, 6))
    ax[0, 0].plot(x_lam1, y_vec1[:, ::-1], lw=0.5)
    ax[0, 0].set_ylim([-0.1, 0.1])
    ax[1, 0].plot(x_lam1, y_vec2[:, ::-1], lw=0.5)
    ax[1, 0].set_ylim([-0.1, 0.1])
    ax[0, 1].plot(x_lam2, y_vec3[:, ::-1], lw=0.5)
    ax[0, 1].set_ylim([-0.1, 0.1])
    ax[1, 1].plot(x_lam2, y_vec4[:, ::-1], lw=0.5)
    ax[1, 1].set_ylim([-0.1, 0.1])

    ax[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0, 1].xaxis.set_major_formatter(plt.NullFormatter())

    ax[0, 0].set_title('PCA vecs with sky, no weights')
    ax[1, 0].set_title('PCA vecs with sky, ivar as weights')
    ax[0, 1].set_title('PCA vecs without sky, no weights')
    ax[1, 1].set_title('PCA vecs without sky, ivar as weights')

    fig.suptitle('First 3 principal vectors (from %s gals)' % number, fontsize=16)
    plt.savefig("plots/sdss/PCA_SDSS_vecs.pdf")
    plt.show()


def calculate_magnitude_r(lam, spec):
    """
    Calcuulate the r-band magnitude of a given spectrum for the DES camera
    :param lam: wavelength grid, log or linear
    :param spec: spectrum
    :return mag_r: r magnitude for the DES camera
    """

    c = 3 * 1e8  # m/s

    # given filter specification, total throughput and a spectrum, calcualte magnitude
    # nu lam = c, nu = c / lam, d_nu = c / lam^2 d_lam
    # TODO: find the right values for sdss or use the filter curves

    l_mean = 6250.  # sdss r-band mean
    delta_l = 1400.  # sdss r-band FWHM
    l_min = l_mean - 0.5 * delta_l
    l_max = l_mean + 0.5 * delta_l

    mask = (lam > l_min) * (lam < l_max)
    lam = lam[mask]
    d_lam = lam[1:] - lam[:-1]  # A
    lam = (lam[1:] + lam[:-1]) / 2  # A
    spec = spec[mask]
    spec = (spec[1:] + spec[:-1]) / 2  # erg/s/cm^2/A
    # https://en.wikipedia.org/wiki/AB_magnitude
    int1 = np.sum(spec * lam / (c * 1e10) * d_lam)
    int2 = np.sum(1. / lam * d_lam)  # 1/s=Hz
    # write this again
    mag_r = -2.5 * np.log10(int1 / int2) - 48.6

    return mag_r


def check_mag_r(lam, spec, mag_r):
    """
    Calculate the r band magnitude and compares it to given magnitude, prints a warning if they differ by more than 5%
    :param lam: wavelength grid
    :param spec: spectrum
    :param mag_r: magnitude to compare with
    :return mag_r_calc: calculated mag_r of the spectrum
    """

    mag_r_calc = calculate_magnitude_r(lam, spec)

    return mag_r_calc