# Copyright 2013 ETHZ.ch Lukas Gamper <lukas.gamper@usystems.ch>
from uspec.utils.utils import Struct

global_ctx = None

def ctx():
    """
    Returns the current global namespace context.
    
    :return: reference to the context module
    
    """
    global global_ctx
    if(global_ctx is None):
        global_ctx = _create_ctx()
        
    return global_ctx

def _create_ctx():
    return Struct()