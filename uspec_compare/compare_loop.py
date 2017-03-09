#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file returns the angles and areas of the 2-sigma ellipses of the coefficient correlation distribution
# or the coefficients themselves
#
# System imports
import numpy as np
import astropy.io.fits as aif
import hope

# Local modules
import plot_coeffs as pcs


def compare(lam, spec, distance='sigma', what='whatever'):
    '''
    Get the data (coefficient scatter angle and area) to use in the control loop

    :param lam: wave length grid of the spectra
    :param spec: galaxy spectra to conduct the pca with
    :param distance: which distance is used in control loop (if default then 2 sigma ellipses are used to get angles and areas)
    :param what: name of the spectra set, here it doesn't matter (only if one does the scatter plot)
    :return aaa: angles and area of the 2-sigma ellipse of the scatter plots from pca coefficients wrt kcorrect
    '''

    kcorrect_template = aif.open('uspec/data/k_nmf_derived.newdefault.fits', memmap=True)

    # kcorrect components and wave length grid
    k_lam = kcorrect_template[11].data
    templates_v0 = np.array(kcorrect_template[5].data, dtype='float32')
    templates = np.zeros((5,len(lam)))
    for idx in range(len(templates)):
        templates[idx] = np.interp(lam, k_lam, templates_v0[idx])

    # conduct pca of uspec spectra with kcorrect components
    coeff = np.zeros((len(spec),5))
    for idx in range(len(spec)):
        coeff[idx] = np.linalg.lstsq(templates.T, spec[idx])[0]

    if (distance=='sigma'):
        # scatter plot of 2 coefficients at a time of fist five sdss PCA coefficients wrt. PCA of uspec (change and to wrt)
        angles, areas = pcs.coeff_scatter(coeff, what, graphic=False)
        aaa = np.array([angles,areas])
        return aaa

    else:
        return coeff


def mmd(x, y, sigma=1e11):
    '''
    Maximum Mean Discrepancy distance

    :param x: first data set. Expected shape (N, ndim)
    :param y: second data set. Expected shape (N, ndim)
    :param sigma: Sigma for gaussian RBF. Rule of thumb: Choose sigma such that ||x-y||**2 / (2 * sigma**2)) equals 1
    for the median distance between x and y
    :returns mmd: the distance between both data sets
    '''

    # Cast to correct dtype
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    sigma = np.float64(sigma)
    lamb = np.float64(-1. / 2. / sigma ** 2)
    n_samples = len(x)

    # Compute contributions to MMD distance
    kxx_sum = _gauss_kernel_vv(x, n_samples, lamb)
    kyy_sum = _gauss_kernel_vv(y, n_samples, lamb)
    kxy_sum = np.zeros(n_samples)
    kyx_sum = np.zeros(n_samples)
    _gauss_kernel_vw(x, y, n_samples, lamb, kxy_sum, kyx_sum)

    kxy_diag = np.exp(lamb * np.sum((x - y) ** 2, axis=1))

    h = kxx_sum + kyy_sum - kxy_sum - kyx_sum + 2 * kxy_diag

    mmd = np.sum(h) / n_samples / (n_samples - 1)
    print(mmd)

    return mmd


@hope.jit
def _gauss_kernel_vw(v, w, n, lamb, out0, out1):
    for i in range(n):
        for j in range(n):
            kvw = np.exp(lamb * np.sum((v[i, :] - w[j, :]) ** 2))
            out0[i] += kvw
            out1[j] += kvw


@hope.jit
def _gauss_kernel_vv(v, n, lamb):
    out = np.zeros(n)
    for i in range(0, n):
        for j in range(i + 1, n):
            kvv = np.exp(lamb * np.sum((v[i, :] - v[j, :]) ** 2))
            out[i] += kvv
            out[j] += kvv
    return out
