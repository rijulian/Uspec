=============================
uspec
=============================

This is an ultra fast spectroscopy simulator that generates realistic observed spectra for a user-specified instrument. The simulator is set up to model data taken in a spectropic survey. 


Features
--------

* Modular structure that allows users to swap in and out different steps in the process of carrying out a spectroscopic survey
* Can be combined with the Ultra Fast Image Simulator (UFIG) 


Prerequisites
-------------

Installing DarkSkySync
--------

http://refreweb.phys.ethz.ch/software/darkskysync/0.1.5/installation.html#how-to-set-up-a-password-less-ssh-login


Installing uspec
----------------

Install the required packages

$ python setup.py develop --user

Now get the data, you can use a different file if you like:

$ cd uspec/data
$ curl -O http://data.sdss3.org/sas/dr10/boss/spectro/redux/v5_5_12/3872/spPlate-3872-55382.fits
$ curl -O http://data.sdss3.org/sas/dr10/boss/spectro/redux/v5_5_12/3872/v5_5_12/spZbest-3872-55382.fits
$ curl -O http://data.sdss3.org/sas/dr10/sdss/spectro/redux/26/0266/spPlate-0266-51602.fits
$ curl -O http://data.sdss3.org/sas/dr10/sdss/spectro/redux/26/0266/spZbest-0266-51602.fits


Running uspec
-------------

Usage:

ivy [options] CONFIG

Examples:

$ ivy uspec.config.sdss
or 
$ ivy uspec.config.bcc
or 
$ ivy uspec.config.bcc_nomr

The possible configurations are:
	bcc: takes BCC coefficients, generates noisy spectra and measures the redshifts
	sdss: takes SDSS coefficients, generates noisy spectra and measures the redshifts
   	bcc_nomr: does not call measure redshift


Useful options
--------------

--debug
	creates additional debug data and plots, will run slower with debug=True. Valid options are: True, false

--sync
	you can disable DarkSkySync using this switch. It's useful if you want to use uspec without internet connections,
	for example on an airplane. The files used by uspec have to be downloaded using DarkSkySync before. sync=False disables
	runs uspec without using DarkSkySync. Valid options are: True, false

--cut
	there are different cuts built into uspec, cmass lowz sdssmain, none. This option can be used to select any of the cuts.
