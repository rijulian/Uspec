#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# according to:
# Copyright (C) 2013 ETH Zurich, lgamper@phys.ethz.ch
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.


# System imports
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np

class CmassFilter(object):
    '''
    Filter that applies the CMASS cut

     :param ModelMag_g: ModelMag g or Mag g
     :param ModelMag_r: ModelMag r or Mag r
     :param ModelMag_i: ModelMag i or Mag i
     :param cModelMag_i: cModelMag i or Mag i
     :param fiber2Mag_i: Fiber2Mag i or Mag i
     :param z_gal: z

     :returns mask: boolean mask to be used on input parameters

    '''   
    def __init__(self, ModelMag_g, ModelMag_r, ModelMag_i, cModelMag_i, fiber2Mag_i, z_gal):
        self.ModelMag_g = ModelMag_g
        self.ModelMag_r = ModelMag_r
        self.ModelMag_i = ModelMag_i
        self.cModelMag_i = cModelMag_i
        self.fiber2Mag_i = fiber2Mag_i
        self.z_gal = z_gal
    
    def get_cmass(self):
        mask = self.z_gal < 0.7
        mask *= self.z_gal > 0.4
        mask *= self.cModelMag_i > 17.5
        mask *= self.cModelMag_i < 19.9
        mask *= self.fiber2Mag_i < 21.5
        mask *= (self.ModelMag_r - self.ModelMag_i) < 2.
        
        d = np.zeros_like(self.ModelMag_r)
        
        d = self.ModelMag_r - self.ModelMag_i - (self.ModelMag_g - self.ModelMag_r) / 8.
        
        mask *= d > 0.55
        mask *= self.cModelMag_i < (19.86 + 1.6 * (d - 0.8))
        return mask