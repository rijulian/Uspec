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

Now get the data

$ cd uspec/data


Running uspec
-------------

Usage:

ivy [options] CONFIG

Examples:

$ ivy uspec.config.ufig_nomr
