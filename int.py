import os
import sys
from time import localtime, strftime
from osgeo import gdal
from osgeo import osr
import numpy as np



def print_time(Action):
        '''Print local time during some given action'''
	Time=strftime("%a, %d %b %Y %H:%M:%S", localtime())
	print(Action+" at "+Time)

	
def choose_prior():
        '''Choose prior field which gives the edges of the file to improve interpolation at the edges of the field'''
	print("prior:")
	
	
def read(FileName, GetProj=False):
        '''Read GeoTiff and save the information of the first band as a numpy array.
        input:  FileName: string containing path to file which should be read
        output: Array: 2d numpy array containing the information of the first band of the GeoTiff
                GeoT:  geographic transformation
                Proj:  projection'''
	print("read file "+filename)
	try:
		DataSet = gdal.Open(FileName, gdal.GA_ReadOnly)
	except IOError as (errno, strerror):
		print "I/O error({0}): {1}".format(errno, strerror)
	assert DataSet.RasterCount == 1, "More than one band in GeoTiff"
	Band = DataSet.GetRasterBand(1)
	Array = Band.ReadAsArray()
        if GetProj=True:
                GeoT = DataSet.GetGeoTransform()
                Proj = osr.SpatialReference()
                Proj.ImportFromWkt(DataSet.GetProjectionRef())    
                return Array, GeoT, Proj
        else:
                return Array


def write(Array, FileName, GeoT, Proj):
        '''Write numpy array data as GeoTiff.
        input: Array:     2d numpy array
               FileName: path to the output file (GeoTiff) ''' 
	print("write file "+FileName)
        Driver = gdal.GetDriverByName('GTiff')
        M = Array.shape[0]
        N = Array.shape[1]
        Type = gdal.GDT_Float32
        DataSet = Driver.Create(FileName, M, N, 1, Type) # the '1' is for band 1.
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Proj.ExportToWkt())
        DataSet.GetRasterBand(1).WriteArray(Data)
        
	
def add_edge(Array, Edge):
    '''Add edges from a prior field to an array
    input: Array: 2d numpy array of the size mxn
           Edge : 2d numpy array of the size mxn
    output: 2d numpy array of the size mxn which contains the edges of the array edge. The rest of the data is the same as in the array array.'''
    N = Array.shape[0] - 1
    M = Array.shape[1] - 1
    Array[0, :] = Edge[0, :]
    Array[:, 0] = Edge[:, 0]
    Array[N, :] = Edge[N, :]
    Array[:, M] = Edge[:, M]
    return Array
    	

def main():
	print_time("Computation started")
	InFile = sys.argv[1]
	OutFile = sys.argv[2]
        DirName=OutFile.rpartition("/")[0]
	#Proj = sys.argv[3]
	PriorFile = choose_prior()
	Input, InputGeoT, InputProj = read(Infile, Proj=True)
	Prior = read(PriorFile)
	Input = add_edge(Input, Prior)
	Isnan = np.isnan(Input)
        write(Isnan,DirName+"/isnan.tif")
        os.system("mkdir "+DirName+"/vectorized")
        os.system("gdal_polygonize.py -q "+DirName+"/isnan.tif"+" -f 'ESRI Shapefile'"+DirName+"/vectorized")
        
	
if __name__ == "__main__":
    main()
