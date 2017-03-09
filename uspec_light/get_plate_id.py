#! /usr/bin/env python
# Copyright (C) 2016 ETH Zurich, rijulian@phys.ethz.ch
#
# This file selects the CMASS galaxies of the SDSS files stored in uspec/data/DR13
#
# System imports
from __future__ import print_function, division, absolute_import, unicode_literals

import astropy.io.fits as aif
import numpy as np

def get_plate_id():
    '''
    Create boolean mask of CMASS cut galaxies to be applied to spZbest files and return the corresponding CMASS galaxy
    redshifts, i-band magnitudes and r-band magnitudes

    :return plate_mask: boolean mask to single out CMASS cut galaxies
    :return z_plate: CMASS galaxy redshifts
    :return mag_plate: CMASS galaxy i-band magnitudes
    :return magr_plate: CMASS galaxy r-band magnitudes
    '''

    # load SDSS SpecPhotoAll and PhotoObjAll catalogs (taken from SDSS Skyserver)
    spec_photo_all = 'uspec/data/SDSS_sample/sdss_cmass_spec_petro.fit'
    photo_obj_all = 'uspec/data/SDSS_sample/sdss_cmass_photo.fit'

    spec_data = aif.open(spec_photo_all)
    photo_data = aif.open(photo_obj_all)

    # get all the CMASS galaxy IDs of the SDSS
    data_spec = spec_data[1].data
    data_photo = photo_data[1].data
    cmass_mask = np.in1d(data_spec['specObjID'], data_photo['specObjID'])
    spec_cmass = data_spec[cmass_mask]
    mag_i = data_spec['petroMag_i'][cmass_mask]
    mag_r = data_spec['petroMag_r'][cmass_mask]


    # load SDSS spZbest catalog
    # original source of the plates is:
    # https://data.sdss.org/sas/dr13/eboss/spectro/redux/v5_9_0/

    #ten random plates
    pl_num = np.zeros((10,2))
    
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

    # compare the galaxy IDs of the CMASS galaxies with the ones in the spZbest files and get the corresponding
    # redshifts and magnitudes
    plate_mask = {}
    z_plate = []
    mag_plate = []
    magr_plate = []

    for idx in range(len(pl_num)):
        sdss_spec_spzbest = 'uspec/data/DR13/spZbest-%s-%s.fits'%(int(pl_num[idx][0]),int(pl_num[idx][1]))
        spZbest_plate = aif.open(sdss_spec_spzbest)
        data_spZbest = spZbest_plate[1].data
        plate_mask[idx] = np.in1d(data_spZbest['specObjID'], spec_cmass['specObjID'])
        plate_mask_mag = np.in1d(spec_cmass['specObjID'], data_spZbest['specObjID'])
        z_plate.append(data_spZbest['Z'][plate_mask[idx]])
        mag_plate.append(mag_i[plate_mask_mag])
        magr_plate.append(mag_r[plate_mask_mag])
    
    return plate_mask, z_plate, mag_plate, magr_plate
    