#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@student.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# This program contains some analysing tools for USPEC spectra


# System imports
from __future__ import print_function, division, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as plt


def plot_spec(survey_lam, spec, vecs, mag_i, z_spec):
    '''
    Plot three spectra with and without noise (incl energy level), the reconstructed spectra and first principal vectors

    :param survey_lam: wave length grid
    :param spec: dictionary with all the different types of spectra
    :param vecs: first principal vectors
    :param mag_i: i-band magnitudes
    :param z_spec: redshifts
    '''

    fig, a = plt.subplots(3, 3, figsize=(20, 6))
    a = a.ravel()

    for idx, ax in enumerate(a):
        if (idx < 3):
            ax.plot(survey_lam, spec['%s' % idx], c='black', lw=0.5)
            ax.plot(survey_lam, spec['%s' % (idx + 9)], 'r-', lw=0.5)
            ax.set_title('Mag_i: %s, z: %s' % (mag_i[idx % 3], z_spec[idx % 3]))
        if (idx > 2):
            ax.plot(survey_lam, spec['%s' % idx], c='black', lw=0.5)
        if (idx > 5):
            ax.set_xlabel('$\lambda \ [\AA]$')
        if (idx == 3):
            ax.set_ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
        if (idx < 6):
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.set_ylim([-2, 5])

    fig.suptitle('Uspec spectra before and after noise and PCA reconstruction', fontsize=16)
    plt.savefig("plots/uspec/spec_uspec_compare.pdf")
    plt.show()

    plt.figure()
    for i in range(len(vecs)):
        plt.plot(survey_lam, vecs[i], lw=0.5, zorder=5 - i, label='Principal vector %s' % (i + 1))
    plt.ylim([-0.1, 0.1])
    plt.legend(loc='best')
    plt.title('Principal vectors')
    plt.savefig("plots/uspec/principal_vecs.pdf")
    plt.show()


def plot_ivar(survey_lam, nine_gals_ivar, mag_i, z_nine_gals):
    '''
    Plot three inverse variances corresponding to three galaxy spectra

    :param survey_lam: wave length grid
    :param nine_gals_ivar: inverse variances
    :param mag_i: i-band magnitudes
    :param z_nine_gals: redshifts
    '''
    fig, a = plt.subplots(1, 3, figsize=(20, 5))
    a = a.ravel()
    for idx, ax in enumerate(a):
        ax.plot(survey_lam, nine_gals_ivar[idx], c='black', lw=0.5)
        ax.text(3600, 60, r'$mag_i$: %s, $z$: %s' % (mag_i[idx], z_nine_gals[idx]))
        ax.set_xlabel('$\lambda \ [\AA]$')
        ax.set_ylim([0, 70])
    fig.suptitle('Inverse Variances', fontsize=16)
    plt.savefig("plots/uspec/multiple_ivar.pdf")
    plt.show()


def ivar_val(survey_lam, spec_pca, spec_nonoise, noise, ivar_pca, num_gals):
    '''
    Compare the observed spectrum of a galaxy with its eigenspectrum, the noise in its spectrum and its inverse variance

    :param survey_lam: wave length grid
    :param spec_pca: noisy galaxy spectra
    :param spec_nonoise: eigenspectra of the galaxies
    :param noise: noise levels
    :param ivar_pca: inverse variances
    :param num_gals: number of galaxies available
    '''
    ran_num = np.random.choice(num_gals)
    fig, ax = plt.subplots(2, 2, figsize=(16, 6))
    ax[0, 0].plot(survey_lam, spec_pca[ran_num], c='black', lw=0.5)
    ax[0, 0].set_ylim([-2, 5])
    ax[0, 1].plot(survey_lam, spec_nonoise[ran_num], c='black', lw=0.5)
    ax[0, 1].set_ylim([-2, 5])
    ax[1, 0].plot(survey_lam, noise[ran_num], c='black', lw=0.5)
    ax[1, 1].plot(survey_lam, ivar_pca[ran_num], c='black', lw=0.5)
    ax[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0, 1].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0, 0].set_title('Observed spectrum')
    ax[0, 1].set_title('Galaxy spectrum wo. Poisson')
    ax[1, 0].set_title('Noise wo. Poisson')
    ax[1, 1].set_title('Ivar')
    ax[1, 0].set_xlabel('$\lambda \ [\AA]$')
    ax[1, 1].set_xlabel('$\lambda \ [\AA]$')
    ax[0, 0].set_ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
    fig.suptitle('Uspec Validation', fontsize=16)
    plt.savefig("plots/uspec/Ivar_Validation.pdf")
    plt.show()

    plt.figure()
    plt.title('Inverse Variance')
    plt.plot(survey_lam, ivar_pca[ran_num], c='black', lw=0.5)
    plt.xlabel('$\lambda \ [\AA]$')
    plt.savefig("plots/uspec/Ivar.pdf")
    plt.show()


def plot_kcorrect(survey_lam, lam, templates):
    '''
    Plot the kcorrect templates interpolated to survey wave length grid

    :param survey_lam: wave length grid of survey
    :param lam: wave length grid of kcorrect templates
    :param templates: kcorrect templates
    '''
    template1 = np.interp(survey_lam, lam, templates[0])
    template2 = np.interp(survey_lam, lam, templates[1])
    template3 = np.interp(survey_lam, lam, templates[2])
    template4 = np.interp(survey_lam, lam, templates[3])
    template5 = np.interp(survey_lam, lam, templates[4])
    plt.figure()
    plt.plot(survey_lam, template1, label=r'Kcorrect principal component 1')
    plt.plot(survey_lam, template2 * 1e-2, label=r'Kcorrect principal component 2 ($\cdot 10^{-2}$)')
    plt.plot(survey_lam, template3, label=r'Kcorrect principal component 3')
    plt.plot(survey_lam, template4, label=r'Kcorrect principal component 4')
    plt.plot(survey_lam, template5, label=r'Kcorrect principal component 5')
    plt.legend(loc='best')
    plt.title('Kcorrect principal components')
    plt.xlabel('$\lambda \ [\AA]$')
    plt.savefig("plots/uspec/Kcorrect.pdf")
    plt.show()


def plot_random_specs(survey_lam, spec_pca, spec_nonoise, noise, z_gals, magnitude_i, num_gals):
    '''
    Plot nine random galaxy spectra with and without noise and the corresponding noise level

    :param survey_lam: wave length grid
    :param spec_pca: noisy galaxy spectra
    :param spec_nonoise: noiseless galaxy spectra
    :param noise: noise levels
    :param z_gals: redshifts
    :param magnitude_i: i-band magnitudes
    :param num_gals: number of available galaxies
    '''
    ran_idx = np.random.choice(num_gals, 9)
    plt.figure(figsize=(25, 10))
    idxs = 1
    for idx in ran_idx:
        plt.subplot(3, 3, idxs)
        plt.title(r'$z$ = %s , $mag_i$ = %s' % (z_gals[idx], magnitude_i[idx]))
        plt.plot(survey_lam, spec_pca[idx], color="k", zorder=2, linewidth=0.5, label="Noisy spectrum of galaxy")
        plt.plot(survey_lam, spec_nonoise[idx], color="r", zorder=3, linewidth=1.,
                 label="Noise-free spectrum of galaxy")
        plt.plot(survey_lam, noise[idx], color="lightskyblue", zorder=1, linewidth=0.5,
                 label="Noise in galaxy spectrum")
        plt.ylim([-2., 10.])
        plt.legend(loc='best', prop={'size': 6})
        plt.ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
        if (idxs > 6):
            plt.xlabel('$\lambda \ [\AA]$')
        idxs += 1
    plt.suptitle('Random Uspec CMASS galaxies', fontsize=18)
    plt.savefig("plots/uspec/spec.pdf")
    plt.show()


def plot_readnoise_sky(survey_lam, read_noise, sky, num_gals):
    '''
    Plot the read noise and observed sky of a random galaxy

    :param survey_lam: wave length grid
    :param read_noise: read noises
    :param sky: skys
    :param num_gals: number of available galaxies
    '''
    ran_idx = np.random.choice(num_gals)
    plt.figure()
    plt.title('Read Noise')
    plt.plot(survey_lam, read_noise[ran_idx], color="k", linewidth=0.5)
    plt.ylim([-100, 100])
    plt.ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
    plt.xlabel('$\lambda \ [\AA]$')
    plt.savefig("plots/readnoise.pdf")
    plt.show()
    plt.figure()
    plt.title('Observed Sky Spectrum')
    plt.plot(survey_lam, sky[ran_idx], color="k", linewidth=0.5)
    plt.ylabel('$Flux \ (10^{-17})[erg/s/cm^2/ \AA]$')
    plt.xlabel('$\lambda \ [\AA]$')
    plt.savefig("plots/uspec/sky.pdf")
    plt.show()