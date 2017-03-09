#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# This file conducts PCAs with a given sample

# System imports
from __future__ import print_function, division, absolute_import, unicode_literals

# External modules
import numpy as np
from wpca import PCA, WPCA
import matplotlib.pyplot as plt


###########
#   PCA   #
###########


def pca_sample(ThisPCA, X, weights=None, ncomp=1):
    '''
    Conduct PCA of a galaxy spectra set

    :param ThisPCA: Type of PCA (eg. PCA, WPCA, EMPCA)
    :param X: galaxy spectra set
    :param weights: weights used in the galaxy spectra PCA (only possible for WPCA), if specified
    :param ncomp: number of principal components

    :return Y.T: reconstructed galaxy spectra
    :return vecs: first ncomp principal vectors
    :return coeffs: first ncomp coeffs of all galaxy spectra
    '''
    # Compute the standard/weighted PCA
    if weights is None:
        kwds = {}
    else:
        kwds = {'weights': weights}

    # Compute the PCA vectors & variance
    pca = ThisPCA(n_components=10).fit(X, **kwds)
    vecs = pca.components_[:ncomp].T

    # Reconstruct the data using the PCA model
    Y = ThisPCA(n_components=ncomp).fit_reconstruct(X, **kwds)

    # Compute the first ncomp coefficients
    coeffs = ThisPCA(n_components=ncomp).fit_transform(X, **kwds)

    return Y.T, vecs, coeffs


def plot_results(ThisPCA, x_val, X, num_gals, what, weights=None, ncomp=1):
    '''
    Conduct PCA of a galaxy spectra set and plot the result

    :param ThisPCA: Type of PCA (eg. PCA, WPCA, EMPCA)
    :param x_val: wave length grid
    :param X: galaxy spectra set
    :param num_gals: number of galaxies to conduct PCA with
    :param weights: weights used in the galaxy spectra PCA (only possible for WPCA), if specified
    :param ncomp: number of principal components
    '''
    # Choose random galaxy to reconstruct
    ran_num = np.random.choice(num_gals)

    # Compute the standard/weighted PCA
    if weights is None:
        kwds = {}
    else:
        kwds = {'weights': weights}

    # Compute the PCA vectors & variance
    pca = ThisPCA(n_components=10).fit(X, **kwds)

    # Reconstruct the data using the PCA model
    Y = ThisPCA(n_components=ncomp).fit_reconstruct(X, **kwds)
    # Create the plots
    fig, ax = plt.subplots(2, 2, figsize=(16, 6))
    ax[0, 0].plot(x_val, X[ran_num].T, c='black', lw=0.5)
    ax[0, 0].set_ylim([-2, 5])
    ax[1, 1].plot(x_val, Y[ran_num].T, c='black', lw=0.5)
    ax[1, 1].set_ylim([-2, 5])

    ax[0, 1].plot(pca.components_[:ncomp].T, lw=0.5)
    ax[0, 1].set_ylim([-0.1, 0.1])

    ax[1, 0].plot(np.arange(1, 11), pca.explained_variance_ratio_)
    ax[1, 0].set_xlim(1, 10)
    ax[1, 0].set_ylim(0, None)
    ax[1, 0].text(3, 0.15, 'Variance ratio covered by the first 5 vectors: %s'%(int(np.sum(pca.explained_variance_ratio_[:5]) * 1000) / 1000.), style='italic')

    ax[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0, 1].xaxis.set_major_formatter(plt.NullFormatter())

    ax[0, 0].set_title('Input Data (random galaxy)')
    ax[0, 1].set_title('First {0} Principal Vectors'.format(ncomp))
    ax[1, 1].set_title('Reconstructed Data ({0} components)'.format(ncomp))
    ax[1, 0].set_title('PCA variance ratio')
    ax[1, 0].set_xlabel('principal vector')
    ax[1, 0].set_ylabel('proportion of total variance')

    fig.suptitle(ThisPCA.__name__ + ' of ' + str(num_gals) + ' CMASS cut galaxies', fontsize=16)
    plt.savefig("%s/PCA_%s.pdf"%(what,what))
    plt.show()

    plt.figure()
    plt.plot(x_val, X[ran_num].T, label='Observed SDSS spectrum', lw=0.5)
    plt.plot(x_val, Y[ran_num].T, label='Reconstructed with PCA', lw=0.5)
    plt.ylim([-2, 5])
    plt.legend(loc='best')
    plt.savefig("%s/PCA_%s_compare.pdf"%(what,what))
    plt.show()

    plt.figure()
    for i in range(ncomp):
        plt.plot(x_val, pca.components_[i].T, lw=0.5, zorder=5-i, label='Principal component %s'%(i + 1))
    plt.legend(loc='best')
    plt.xlabel('$\lambda \ [\AA]$')
    plt.title('Principal vectors')
    plt.savefig("%s/PC_%s.pdf"%(what,what))
    plt.show()
