# Copyright (c) 2016 ETH Zurich, Institute of Astronomy, Cosmology Research Group

"""
@author: Joerg Herbel
"""

# Imports
from __future__ import print_function, division, absolute_import, unicode_literals
import ivy


def run_uspec_from_config(config, **kwargs):
    """
    Run USPEC using an (importable) configuration file. Additional parameters can be provided as keyword arguments.
    :param config: Importable configuration file.
    :param kwargs: Additional parameter values.
    :return: Ivy-context created by USPEC
    """
    config = dict(ivy.loadConfigs(config))
    config.update(**kwargs)
    ctx = run_uspec_from_dict(config)
    return ctx


def run_uspec_from_dict(dictionary):
    """
    Run USPEC using a dictionary as configuration.
    :param dictionary: Configuration dictionary.
    :return: Ivy-context created by USPEC
    """
    ctx = ivy.execute(dictionary)
    return ctx