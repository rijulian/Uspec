# Copyright (c) 2013 ETH Zurich, Institute for Astronomy

'''
Created on Feb 12, 2014

@author: jakeret
'''

from __future__ import print_function, division, absolute_import, unicode_literals

import UserDict

# In Python 2.7 still, `DictMixin` is an old-style class; thus, we need
# to make `Struct` inherit from `object` otherwise we loose properties
# when setting/pickling/unpickling

class Struct(object, UserDict.DictMixin):
    """
    A `dict`-like object, whose keys can be accessed with the usual
    '[...]' lookup syntax, or with the '.' get attribute syntax.

    Examples::

      >>> a = Struct()
      >>> a['x'] = 1
      >>> a.x
      1
      >>> a.y = 2
      >>> a['y']
      2

    Values can also be initially set by specifying them as keyword
    arguments to the constructor::

      >>> a = Struct(z=3)
      >>> a['z']
      3
      >>> a.z
      3

    Like `dict` instances, `Struct`s have a `copy` method to get a
    shallow copy of the instance:

      >>> b = a.copy()
      >>> b.z
      3

    """
    def __init__(self, initializer=None, **extra_args):
        if initializer is not None:
            try:
                # initializer is `dict`-like?
                for name, value in initializer.items():
                    self[name] = value
            except AttributeError:
                # initializer is a sequence of (name,value) pairs?
                for name, value in initializer:
                    self[name] = value
        for name, value in extra_args.items():
            self[name] = value

    def copy(self):
        """Return a (shallow) copy of this `Struct` instance."""
        return Struct(self)

    # the `DictMixin` class defines all std `dict` methods, provided
    # that `__getitem__`, `__setitem__` and `keys` are defined.
    def __setitem__(self, name, val):
        self.__dict__[name] = val
    def __getitem__(self, name):
        return self.__dict__[name]
    def keys(self):
        return self.__dict__.keys()