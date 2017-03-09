#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@student.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.


# System imports
from __future__ import print_function, division, absolute_import, unicode_literals
import numpy as np




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

        # store spectra and variables in arrays and define other variables
        num = len(ctxList[0].spec_observed)
        quarters_len = int(len(par.survey_lam)-len(par.survey_lam)/4)
        spec_pca = np.zeros((len(ctxList),num))

        for idx in range(len(ctxList)):
            spec_pca[idx] = ctxList[idx].spec_observed

        # save in context for control loop to use
        par.spectra = spec_pca[:,:quarters_len]
        par.lam = par.survey_lam[:quarters_len]


    def __str__(self):
        return "reduce and plot plugin"
