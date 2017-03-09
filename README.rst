=============================
Uspec
=============================

This is the control loop for the ultra fast spectroscopy simulator that generates realistic observed spectra for a user-specified instrument. The simulator is set up to model data taken in a spectropic survey. 



Installing uspec
----------------

Install the required packages

$ python setup.py develop --user


Running uspec
-------------

Usage:

ivy [options] CONFIG

Examples:

$ ivy uspec.config.ufig_nomr

Control loop:

$ python control_loop
