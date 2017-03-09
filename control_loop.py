#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@student.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# System imports
from __future__ import print_function, division, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as aif
import abcpmc
import time
import seaborn as sns

# Local modules
import uspec_compare.compare_loop as cl
import uspec_light.pca_sdss as ps
import run_uspec as ru


# Run parameters
runs = 10                                   # number of runs to generate particle pools
alpha = 85                                  # alpha-percentile threshold decrease
eps_start = 0.30                            # starting boundry
eps = abcpmc.ConstEps(runs, eps_start)      # constant threshold, alpha-percentile threshold decrease


def get_data():
    '''
    Get the real data set in order to compare the simulated one with it

    :return data: if distance='sigma', returns angles and areas of 2-sigma ellipses, else returns the coefficients themselves
    '''
    lam, spec = ps.sdss_data()
    data = cl.compare(lam, spec, distance='sigma')
    return data


# Actual data to compare USPEC with
data_0 = get_data()

# Theory parameters
theta_0 = np.array([0.75,1.5,0.75])    # best guess for rn (read noise constant), tr (transmission constant), st (sky transmission constant)
bounds = np.array([0.6,0.6,0.6])  # upper and lower bound estimates on the theory parameters


def create_new_sample(theta):
    '''
    Creates new particle to be used in abcpmc, the data set to be examined is the sample's fitting coefficients
    wrt. kcorrect principal vectors

    :param theta: random theory parameters chosen by prior
    :return aaa: if distance='sigma', returns angles and areas of 2-sigma ellipses, else returns the coefficients themselves
    '''

    ctx = ru.run_uspec_from_config("uspec.config.ufig_nomr",theta=theta)

    lam = ctx.params.lam
    spec = ctx.params.spectra

    aaa = cl.compare(lam, spec, distance='sigma')

    return aaa


def dist_measure(aaa, real_aaa):
    '''
    Generic distance measure of the relative angles residual and the relative areas residual

    :param aaa: angles and areas of the simulated data
    :param real_aaa: angles and areas of the real data
    :return dist: the distance between both data sets
    '''

    dist = (np.mean(np.abs((aaa[0]-real_aaa[0])/real_aaa[0]))+np.mean(np.abs((aaa[1]-real_aaa[1])/real_aaa[1])))/2.
    print('distance:',dist)

    return dist



class flat_prior(object):
    '''
    Flat prior

    :param theta_0: Initial best guess parameters
    :param bounds: Best guess bounds on the parameters
    '''
    def __init__(self, theta_0, bounds):
        self.theta = theta_0
        self.bounds = bounds

    def __call__(self, theta=None):
        if theta is None:
            np.random.seed()
            random_theta = np.random.uniform(theta_0-bounds, theta_0+bounds)
            return random_theta
        else:
            if (np.bool(np.prod(np.array(theta > theta_0-bounds)))) and (np.bool(np.prod(np.array(theta < theta_0+bounds)))):
                return 1./np.prod(2.*bounds)
            else:
                return 0.


prior = flat_prior(theta_0, bounds)
sampler = abcpmc.Sampler(N=100, Y=data_0, postfn=create_new_sample, dist=dist_measure, threads=6)


def launch():
    '''
    Launches the abcpmc routine

    :return pool: pool with all the important data (eg. theta parameters)
    '''

    eps = abcpmc.ConstEps(runs, eps_start)

    pools = []
    for pool in sampler.sample(prior, eps):
        print("T: {0}, eps: {1:>.4f}, ratio: {2:>.4f}".format(pool.t, eps(pool.eps), pool.ratio))

        for i, (mean, std) in enumerate(zip(*abcpmc.weighted_avg_and_std(pool.thetas, pool.ws, axis=0))):
            print(u"    theta[{0}]: {1:>.4f} \u00B1 {2:>.4f}".format(i, mean,std))

        eps.eps = np.percentile(pool.dists, alpha) # reduce eps value
        pools.append(pool)

    sampler.close()
    return pools


# Start the abcpmc control loop
t0 = time.time()
pools = launch()
print("took", (time.time() - t0))



######################
#       plots        #
######################


plt.figure()
moments = []
for i in range(len(theta_0)):
    moments = np.array([abcpmc.weighted_avg_and_std(pool.thetas[:,i], pool.ws, axis=0) for pool in pools])
    plt.errorbar(range(runs), moments[:, 0], moments[:, 1],label='%s parameter'%i)
plt.xlim([-.5, runs])
plt.savefig('values.pdf')

import corner
plt.figure()
samples = np.vstack([pool.thetas for pool in pools])
corner.corner(samples)
plt.savefig('values2.pdf')

sns.set_style("white")
np.random.seed()

plt.figure()
distances = np.array([pool.dists for pool in pools]).flatten()
sns.distplot(distances, axlabel="distance")
plt.savefig('distances.pdf')

plt.figure()
eps_values = np.array([pool.eps for pool in pools])
plt.plot(eps_values, label=r"$\epsilon$ values")
plt.xlabel("Iteration")
plt.ylabel(r"$\epsilon$")
plt.legend(loc="best")
plt.savefig('epsilon.pdf')

plt.figure()
acc_ratios = np.array([pool.ratio for pool in pools])
plt.plot(acc_ratios, label="Acceptance ratio")
plt.ylim([0, 1])
plt.xlabel("Iteration")
plt.ylabel("Acceptance ratio")
plt.legend(loc="best")
plt.savefig('alpha.pdf')

hdu_primary = aif.PrimaryHDU()
cols = aif.ColDefs([
    aif.Column(name=b'params', format=str(len(moments[0])) + 'E', array=moments)])
hdu_data = aif.BinTableHDU.from_columns(cols)
hdulist = aif.HDUList([hdu_primary, hdu_data])
hdulist.writeto("params.fits", overwrite=True)

