 # Copyright (c) 2013 ETH Zurich, Institute of Astronomy, Lukas Gamper <lukas.gamper@usystems.ch>

'''
Main executable to run ufig.

Created on Oct 7, 2013
@author: L. Gamper

'''
import uspec, sys, getopt

def run():
	"""
	Called by the entry point script. Delegating call to main()
	"""
	main(*sys.argv[1:])

def main(*argv):
	if(argv is None or len(argv)<1):
		usage()
		return
		
	if argv[-1].find('.') > -1:
		namespace = getattr(__import__('.'.join(argv[-1].split('.')[:-1]), globals(), locals(), [argv[-1].split('.')[-1]], -1), argv[-1].split('.')[-1])
	else:
		namespace = __import__(argv[-1], globals(), locals(), [], -1)

	args = {}
	args.update((name, getattr(namespace, name)) for name in dir(namespace) if not name.startswith("__"))

	cast = {
		'bool': lambda x: bool(x),
		'int': lambda x: int(x),
		'long': lambda x: long(x),
		'float': lambda x: float(x),
		'str': lambda x: x,
		'list': lambda x: x.split(',')
	}

	# overwrite parameters by command line options
	optlist, positional = getopt.getopt(argv, '', [name.replace('_', '-') + '=' for name in args.keys()])
	if len(positional) != 1:
		raise IOError('only one config file is allowed')
	for opt in optlist:
		if opt[0][:2] != '--':
			raise IOError('invalid option name: {:}'.format(opt[0]))
		elif not opt[0][2:].replace('-', '_') in args:
			raise IOError('unknown option: {:}'.format(opt[0][2:]))
		else:
			args[opt[0][2:].replace('-', '_')] = cast[type(args[opt[0][2:].replace('-', '_')]).__name__](opt[1])

	uspec.generate(**args)

def usage():
	"""
	Return usage of the main uspec call and an example.
	"""

	usage = """
	**Ultra fast spectrum generator**
	Copyright (c) 2014 ETH Zurich, Institute of Astronomy
	
	Usage:
	uspec <arguments> <configurations>
	
	example:
	- uspec uspec.config.bcc
	"""
	print usage

if __name__ == "__main__":
	main(*sys.argv[1:])
	