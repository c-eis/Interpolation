import os
import sys
from time import localtime, strftime
from osgeo import gdal
from osgeo import osr
import numpy as np


def print_time(Action):
    """
    Print local time during some given action
    
    Parameters
    ----------
    Action  : str
        ongoing action at printed time
    """
    Time=strftime("%a, %d %b %Y %H:%M:%S", localtime())
    print(Action+" at "+Time)

	
def choose_prior():
    """
    Choose prior field which gives the edges of the file to improve 
    interpolation at the edges of the field"""
    print("prior:")
	
	
def read(FileName, GetProj=False):
    """
    Read GeoTiff and save the information of the first band as a numpy array.
    
    Parameters
    ----------
    FileName  :  str 
        path to file which should be read
        
    Returns
    -------
    Array  :  2d numpy array 
        information of the first band of the GeoTiff
    GeoT   :  str
        geographic transformation
    Proj   :  str
        geographic projection
    """
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
    """
    Write numpy array data as GeoTiff.
    
    Parameter
    ---------
    Array     :  2d numpy array
        array which should be written in GeoTiff
    FileName  :  str
        path of the output file (GeoTiff)
    GeoT      :  str
        geographic transformation
    Proj      :  str
        geographic projection
    """
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
    """
    Add edges from a prior field to an array
    
    Parameter
    ---------
    Array  :  2d numpy array of the size mxn
    Edge   :  2d numpy array of the size mxn
    
    Returns
    -------
    2d numpy array of the size mxn which contains the edges of the array edge. 
    The rest of the data is the same as in the array array.
    """
    N = Array.shape[0] - 1
    M = Array.shape[1] - 1
    Array[0, :] = Edge[0, :]
    Array[:, 0] = Edge[:, 0]
    Array[N, :] = Edge[N, :]
    Array[:, M] = Edge[:, M]
    return Array
    
    
      
def classification(OldFileName,NewFileName):
    with fiona.open(OldFileName, 'r') as Gaps:
        NumPolygons = len(list(Gaps))
        Properties = np.empty([NumPolygons, 3]) # Characteristics of gaps
        Rad = np.empty([NumPolygons, 1]) # Kriging radius
        i = 0
        for Polygon in Gaps:
            P = shape(Polygon['geometry'])
            Properties[i] = properties(P)
            Rad[i] = diameter(P)+500
            i = i + 1
    Properties = np.divide(Properties, np.absolute(Properties).max(0)) 
        # scaling property values
    NumClusters = 150
    if NumClusters > NumPolygons:
        NumClusters = NumPolygons
    Cluster = KMeans(n_clusters=NumClusters)
    Cluster.fit(Properties)
    Classes = Cluster.predict(Properties)
    for j in xrange(NumClusters):
        Rad[Classes == j] = np.max(Rad[Classes == j])
    i = 0
    with fiona.open(OldFileName, 'r') as Gaps:
        Schema = Gaps.schema
        Schema['properties'].update({u'krig': 'int'})
        with fiona.open(NewFileName, 'w', 'ESRI Shapefile', Schema, Gaps.crs)\
                as NewGaps:
            for Polygon in Gaps:
                Polygon['properties'].update({u'KrigingRadius': int(Rad[i])})
                i = i + 1
                NewGaps.write({'properties': Polygon['properties'],
                           'geometry': Polygon['geometry']})


def extent(FileName):
    Data = gdal.Open(FileName, GA_ReadOnly)
    GeoT = Data.GetGeoTransform()
    MinX = GeoT[0]
    MaxY = GeoT[3]
    MaxX = MinX + GeoT[1] * data.RasterXSize
    MinY = MaxY + GeoT[5] * data.RasterYSize
    Data = None
    return [MinX, MinY, MaxX, MaxY]
    
    
def main():
    print_time("Computation started")
    InFile = sys.argv[1]
    OutFile = sys.argv[2]
    DirName=OutFile.rpartition("/")[0]
    PriorFile = choose_prior()
    Input, InputGeoT, InputProj = read(Infile, Proj=True)
    Prior = read(PriorFile)
    print_time("Add edges")
    Input = add_edge(Input, Prior)
    print_time("Vectorization")
    Isnan = np.isnan(Input)
    IsnanFile = DirName+"/isnan.tif"
    write(Isnan,IsnanFile)
    VecFile = DirName+"/vectorized"
    os.system("mkdir "+VecFile)
    os.system("gdal_polygonize.py -q "+IsnanFile+" -f 'ESRI Shapefile'"+VecFile)
    print_time("Computation of kriging radii")
    RadiusSHPFile = DirName+"/kriging_radius"
    RadiusFile = DirName+"/kriging_radius.tif"
    classification(VecFile, RadiusSHPFile)
    Extent = extent(InFile)
    Extent = str(Extent[0], Extent[1], Extent[2], Extent[3])
    os.system("gdal_rasterize -q -a krig -l "+RadiusSHPFile+" -a_srs "+InputProj+" -tr 250 250 -te "+Extent+" "+RadiusSHPFile+"/kriging_radius.shp "+RadiusFile)
    print_time("Start Kriging")
    #os.system("./kriging.R $surface_gaps_filled.nc  $kriging_radius_inside.nc $kriging_values_inside.nc $interim_inside.nc $kriging_variance_inside.nc $dir")
    
	
if __name__ == "__main__":
    main()
