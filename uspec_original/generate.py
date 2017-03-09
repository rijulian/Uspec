# Copyright (c) 2013 ETH Zurich, Institute of Astronomy, Lukas Gamper <lukas.gamper@usystems.ch>

import time
from uspec.context import ctx

class parameterType(dict):
	def __init__(self, *args, **kwargs):
		dict.__init__(self, *args, **kwargs)
		self.__dict__ = self

# TODO: make test where cat. intimage and psf are run and then do calculation exact and compare!

def generate(**args):

	ctx().parameters = parameterType(args)
	ctx().timings = {}

	ctx().pluginNames = [name.strip() for name in ctx().parameters.plugins]
	ctx().plugins = []
	for name in ctx().pluginNames:
		if name.find('.') > -1:
			module = __import__('.'.join(name.split('.')[:-1]), globals(), locals(), [name.split('.')[-1]], -1)
			ctx().plugins.append(getattr(module, name.split('.')[-1]).plugin())
		else:
			ctx().plugins.append(__import__(name, globals(), locals(), [], -1).plugin())

	# run plugins
	for plugin in ctx().plugins:
		start = time.time()
		plugin()
		ctx().timings[str(plugin)] = time.time() - start
