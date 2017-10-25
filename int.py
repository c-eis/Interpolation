import os
import sys
from time import localtime, strftime
from osgeo import gdal
import numpy as np



def print_time(action):
	t=strftime("%a, %d %b %Y %H:%M:%S", localtime())
	print(action+" at "+t)

	
def choose_prior():
	print("prior:")
	
	
def read(filename):
	print("read file "+filename)
	try:
		dataset = gdal.Open(filename, gdal.GA_ReadOnly)
	except IOError as (errno, strerror):
		print "I/O error({0}): {1}".format(errno, strerror)
	
	band = dataset.GetRasterBand(1)
	array = band.ReadAsArray()
	assert dataset.RasterCount == 1, "More than one band in GeoTiff"
	return array
	
	
def add_edge(array, edge):
    '''Add edges from a prior field to an array'''
    n = array.shape[0] - 1
    m = array.shape[1] - 1
    array[0, :] = edge[0, :]
    array[:, 0] = edge[:, 0]
    array[n, :] = edge[n, :]
    array[:, m] = edge[:, m]
	return array
    	
	
def main():
	print_time("Computation started")
	infile = sys.argv[1]
	outfile = sys.argv[2]
	proj = sys.argv[3]
	priorfile = choose_prior()
	input = read(infile)
	prior = read(prior)
	input = add_edge(input, prior)
	isnan=np.isnan(array)
	
if __name__ == "__main__":
    main()