#Copyright (C) 2016 ETH Zurich, Institute for Astronomy

'''
author: jakeret
recreated: rijulian
'''

from __future__ import print_function, division, absolute_import, unicode_literals


class Plugin(object):
    '''
    classdocs
    '''

    def __init__(self, ctx):
        '''
        Constructor
        '''
        self.ctx = ctx
    
    def getWorkload(self):
        '''
        Divides the workload into pieces to be processed using ivy
    
        '''  
        par = self.ctx.params
        self.ctx.params.numberOfSpectra = len(self.ctx.gal_z)
        
        values = [i for i in range(self.ctx.params.numberOfSpectra)]

        print("spectra count", self.ctx.params.numberOfSpectra)
        
        for value in values:
            ctx = self.ctx.copy()
            ctx.z = self.ctx.gal_z[value]
            del ctx.gal_z
            ctx.mag_r = self.ctx.gal_mag_r[value]
            del ctx.gal_mag_r
            ctx.mag_i = self.ctx.gal_mag_i[value]
            del ctx.gal_mag_i
            ctx.coeffs = self.ctx.gal_coeffs[value]
            del ctx.gal_coeffs
            ctx.galaxy_index = value
            
            yield ctx


    def __str__(self):
        return "map plugin"
