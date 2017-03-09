import numpy as np
from scipy import integrate

def binspec(x_old, y_old, x_new, y_new):
	
	a1 = np.min(x_old)
	a2 = np.max(x_old)
	b1 = np.min(x_new)
	b2 = np.max(x_new)
	
	
	if (b1>=a1 and b2<=a2):
		
		bin_centers = x_new

		Nnew = np.size(bin_centers)
		
		bin_edges = np.ndarray((Nnew+1,), dtype='float32')
		bin_edges[1:Nnew] = (bin_centers[1:]+bin_centers[:-1])*.5
		bin_edges[0] = max(x_old[0],bin_centers[0]-(bin_centers[1]-bin_centers[0])*.5)
		bin_edges[-1] = min(x_old[-1],bin_centers[-1]+(bin_centers[-1]-bin_centers[-2])*.5)

		bin_widths = bin_edges[1:] - bin_edges[:-1]

		y_old_interp = np.interp(bin_edges,x_old,y_old)
		
		insert_index = []
		for j in range(Nnew+1):
			i = np.sum(x_old < bin_edges[j])
			insert_index.append(i)
		x_old_refined = np.insert(x_old,insert_index,bin_edges)
		y_old_refined = np.insert(y_old,insert_index,y_old_interp.astype('float32'))

		for j in range(Nnew):
			istart = insert_index[j] + j
			istop = insert_index[j+1] + (j+1)
			x_old_slice = x_old_refined[range(istart,istop+1)]
			y_old_slice = y_old_refined[range(istart,istop+1)]
			y_new[j] = integrate.trapz(y_old_slice,x_old_slice)
		return True
	else:
		y_new = []
		print "Can't rebin, redshift is too large"
		return False

def smoothspec(loglam,flux,vdisp):

	nlambda		  = len(loglam)
	pixsize		  = np.abs(np.log(10.)* 2.99792e5* (loglam[nlambda-1]-loglam[0])/float(nlambda))
	smoothing	  = float(vdisp)/float(pixsize)								 # pixels
	npix		  = long(4.0* np.ceil(smoothing))* long(2)+3
	klam		  = np.arange(float(npix))- float(npix- 1.)/ 2.
	kernel		  = np.exp(-0.5* (klam/ smoothing)**2)/np.sqrt(2.*np.pi)/smoothing
	kernel		  = kernel/np.sum(kernel)
	smoothed_spec = np.convolve(flux,kernel,mode='same')#,/edge_truncate)	python code can't do this yet.# "valid" mode gives too short

	return smoothed_spec