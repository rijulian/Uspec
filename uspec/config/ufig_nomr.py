'''
Created in 2016
@author: Julian Riebartsch

Define workflow and input parameters to make a spectrum with ufig
catalogs as input.
'''

from ivy.plugin.parallel_plugin_collection import ParallelPluginCollection

# ==================================================================
# P L U G I N S
# ==================================================================

backend = "sequential"
#backend = "multiprocessing"
#cpu_count = 1

# Plugins to render an image
plugins = [
           "uspec.plugins.import_ufig_cat",
           ParallelPluginCollection([
                                    "uspec.plugins.redshift_spec",
                                    "uspec.plugins.obsSpecUfig"],
                                    "uspec.plugins.spectra_map_ufig",
                                    "uspec.plugins.spectra_ufig"),
           "ivy.plugin.show_stats"
]

survey = 'ufig'
