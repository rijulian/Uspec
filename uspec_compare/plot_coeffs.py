#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file determines angles and areas of the scatter plot of pca coefficients
#
# System imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def coeff_scatter(coeffs, name, graphic=True ,comp=False, coeffs2=None, name2=None, kcorrect=False):
    '''
    Scatter plot of 2 coefficients at a time of first five PCA coefficients of a set and potentially of a second one
    and the corresponding 2 sigma confidence ellipses

    :param coeffs: first set of first five PCA coefficients
    :param name: name of the first set
    :param graphic: if 'False' this code only returns the angles and areas of the 2 sigma ellipse ('True' by default)
    :param comp: if 'True' compare first set of coefficients to the second one
    :param coeffs2: second set of first five PCA coefficients
    :param name2: name of the second set
    :param kcorrect: if "True" both coefficients are done wrt kcorrect pricipal vectors and the angles and areas of the ellipses are returned
    :return angles: angles (degrees) between the lower and higher coefficients (if theta=0, the long side is the lower coeff side)
    :return area: area of the ellipses
    '''
    angles = np.zeros(10)
    area = np.zeros(10)
    if (graphic==True):
        fig, a = plt.subplots(3, 4, figsize=(15, 15),sharex=True, sharey=True)
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]
    if (graphic == True):
        a = a.ravel()
    i = 0
    j = 1
    if (graphic == True):
        for idx, ax in enumerate(a):
            if (idx > 9):
                ax.axis('off')
            else:
                ax.scatter(coeffs[:,i], coeffs[:,i+j],marker=".", alpha=0.25, lw=0, color='maroon')
                cov = np.cov(coeffs[:, i], coeffs[:, i + j])
                vals, vecs = eigsorted(cov)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
                w, h = 2. * 2. * np.sqrt(vals)
                ell = Ellipse(xy=(np.mean(coeffs[:,i]), np.mean(coeffs[:,i+j])), width=w, height=h, angle=theta)
                ell.set_facecolor('none')
                ell.set_edgecolor('maroon')
                ax.add_artist(ell)
                if (comp == True):
                    ax.scatter(coeffs2[:, i], coeffs2[:, i + j], marker=".", alpha=0.25, lw=0, color='royalblue')
                    cov2 = np.cov(coeffs2[:, i], coeffs2[:, i + j])
                    vals2, vecs2 = eigsorted(cov2)
                    theta2 = np.degrees(np.arctan2(*vecs2[:, 0][::-1]))
                    w2, h2 = 2. * 2. * np.sqrt(vals2)
                    ell2 = Ellipse(xy=(np.mean(coeffs2[:, i]), np.mean(coeffs2[:, i + j])), width=w2, height=h2,
                                   angle=theta2)
                    ell2.set_facecolor('none')
                    ell2.set_edgecolor('royalblue')
                    ax.add_artist(ell2)
                ax.set(adjustable='box-forced', aspect='equal')
                ax.set_xlabel('%s. Coefficient'%(i+1))
                ax.set_xlim([-50,50])
                ax.set_ylim([-50,50])
                ax.set_ylabel('%s. Coefficient'%(i+j+1))
                if (idx < 8):
                    ax.xaxis.set_major_formatter(plt.NullFormatter())
                if (idx == 3 or idx == 6 or idx == 8):
                    i += 1
                    j = 1
                else:
                    j += 1
        if (comp == True):
            plt.suptitle(r'Coefficient scatter plot (%s (red) and %s (blue)) and 2$\sigma$ confidence ellipses' %(name, name2), fontsize=25)
            if (kcorrect == True):
                plt.savefig("../plots/scatter_coeffs_%s_%s_%s.pdf"%(name, name2, 'kcorrect'))
                plt.show()
            else:
                plt.savefig("plots/plots_both/scatter_coeffs_%s_%s.pdf"%(name, name2))
                plt.show()
        else:
            plt.suptitle(r'Coefficient scatter plot (%s) and 2$\sigma$ confidence ellipses'%name, fontsize=25)
            plt.savefig("plots/scatter_coeffs_%s.pdf"%name)
            plt.show()
            plt.close()
    else:
        for idx in range(10):
            cov = np.cov(coeffs[:, i], coeffs[:, i + j])
            vals, vecs = eigsorted(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            angles[idx] = theta
            w, h = 2. * 2. * np.sqrt(vals)
            area[idx] = np.pi * w * h / 4.
            if (idx == 3 or idx == 6 or idx == 8):
                i += 1
                j = 1
            else:
                j += 1
        return angles, area
