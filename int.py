from osgeo import gdal
import sys
from time import localtime, strftime
from osgeo import osr
import numpy as np
import os
import fiona
from shapely.geometry import shape, LineString
from sklearn.cluster import KMeans
import subprocess
from subprocess import check_output
import platform



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

    
def choose_prior(FileName, EPSG, Platform):
    """
    Choose prior field which gives the edges of the file to improve 
    interpolation at the edges of the field
    
    Parameter
    ---------
    
    FileName  : str
        Path of the input file
    Proj      : str
        geographic projection of the input file
        
    Returns
    -------
    
    str
        path of the prior field
    """
    if Platform == "Windows":
        # TODO: correct files in directory, there are only dummy copies yet
        if "vx" in FileName:
            PriorNorth = "C:\\Users\\Tine\\Documents\\AWI\\github\\Interpolation\\prior\\prior_vx.tif"
            PriorSouth = "C:\\Users\\Tine\\Documents\\AWI\\github\\Interpolation\\prior\\rignot_vx.tif"
        elif "vy" in FileName:
            PriorNorth = "C:\\Users\\Tine\\Documents\\AWI\\github\\Interpolation\\prior\\prior_vy.tif"
            PriorSouth = "C:\\Users\\Tine\\Documents\\AWI\\github\\Interpolation\\prior\\rignot_vy.tif"
        else:
            # TODO assertion
            print("FileName must include vx or vy!")
            PriorNorth = None
            PriorSouth = None
    elif Platform == "Mac":
        if "vx" in FileName:
            PriorNorth = "/Users/cluettig/Documents/Daten/prior/prior_vx.tif"
            PriorSouth = "/Users/cluettig/Documents/Daten/prior/rignot_vx.tif"
        elif "vy" in FileName:
            PriorNorth = "/Users/cluettig/Documents/Daten/prior/prior_vy.tif"
            PriorSouth = "/Users/cluettig/Documents/Daten/prior/rignot_vy.tif"
        else:
            # TODO assertion
            print("FileName must include vx or vy!")
            PriorNorth = None
            PriorSouth = None
    elif Platform == "Ollie":
        # TODO: Add Ollie Pathes!
        if "vx" in FileName:
            PriorNorth = "/work/ollie/cluettig/prior/prior_vx.tif"
            PriorSouth = "/work/ollie/cluettig/prior/rignot_vx.tif"
        elif "vy" in FileName:
            PriorNorth = "/work/ollie/cluettig/prior/prior_vy.tif"
            PriorSouth = "/work/ollie/cluettig/prior/rignot_vy.tif"
        else:
            # TODO assertion
            print("FileName must include vx or vy!")
            PriorNorth = None
            PriorSouth = None
    else:
        # TODO: assert
        print("Unknown Prior Pathes! Please add pathes for this platform in choose_prior(FileName, EPSG, Platform).")
        
    DataSet = gdal.Open(PriorNorth, gdal.GA_ReadOnly)
    SRS = osr.SpatialReference()
    SRS.ImportFromWkt(DataSet.GetProjectionRef())
    EPSGPrior=SRS.GetAuthorityCode(None)
    print(EPSG)
    #TODO: other EPSG codes possible
    if EPSGPrior == EPSG:
        PriorFile = PriorNorth
    else:
        PriorFile = PriorSouth
    return PriorFile
    
    
def read(FileName, GetGeoT=False, GetWKT=False, GetEPSG=False, GetProj4=False):
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
    print("read file "+FileName)
    #try:
    DataSet = gdal.Open(FileName, gdal.GA_ReadOnly)
    #except IOError as (errno, strerror):
    #    print "I/O error({0}): {1}".format(errno, strerror)
    assert DataSet.RasterCount == 1, "More than one band in GeoTiff"
    Band = DataSet.GetRasterBand(1)
    Array = Band.ReadAsArray()
    Output = [Array, None, None, None, None]
    if GetGeoT == True:
        GeoT = DataSet.GetGeoTransform()
        Output[1] = GeoT
    if GetWKT == True or GetEPSG == True or GetProj4 == True:
        SRS = osr.SpatialReference()
        SRS.ImportFromWkt(DataSet.GetProjectionRef())
        Output[2] = SRS
    if GetEPSG == True or GetProj4 == True:
        EPSG=SRS.GetAuthorityCode(None)
        Output[3] = EPSG
    if GetProj4 == True:
        Proj4 = SRS.ExportToProj4()
        Output[4] = Proj4
    return Output


def write(Array, FileName, GeoT, WKT):
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
    WKT      :  str
        geographic projection
    """
    print("write file "+FileName)
    Driver = gdal.GetDriverByName('GTiff')
    N = Array.shape[0]
    M = Array.shape[1]
    Array = Array.astype(int)
    Type = gdal.GDT_Int32
    DataSet = Driver.Create(FileName, M, N, 1, Type) # the '1' is for band 1.
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(WKT.ExportToWkt())
    DataSet.GetRasterBand(1).WriteArray(Array)
        
    
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
  
  
def combine(InputArray, PriorArray):
    """
    Add points of PriorArray to InpuArray if they are NaN there.
    
    Parameter
    ---------
    InputArray  : 2d numpy array of the size mxn
    PriorArray  : 2d numpy array of the size mxn
    
    Returns
    -------
    2d numpy array of the size mxn which contains data from InputArray and data 
    from PriorArray at noData positions.
    """
    InputArray[np.isnan(InputArray)]=PriorArray[np.isnan(InputArray)]
      
    
# TODO: is this reasonable?:
def minimal_diameter(Polygon):
    ''' 
    Computes maximal extend of a polygon in x and y direction and returns the 
    smaller one.
    
    Parameter
    ---------
    
    Polygon:  shapely polygon
        A closed polygon 
        
    Returns
    -------
    
    scalar
        minimal diameter of the polygon  
    '''
    [MinX, MinY, MaxX, MaxY] = Polygon.bounds
    return min(MaxX - MinX, MaxY - MinY)

    
# TODO: is this reasonable?:
def maximal_diameter(Polygon):
    ''' 
    Computes maximal extend of a polygon in x and y direction and returns the 
    greater one.
    
    Parameter
    ---------
    
    Polygon:  shapely polygon
        A closed polygon 
        
    Returns
    -------
    
    scalar
        rough diameter of the polygon  
    '''
    [MinX, MinY, MaxX, MaxY] = Polygon.bounds
    Line = LineString([(MinX, MinY), (MaxX, MaxY)])
    return Line.length

 
# TODO: is this reasonable?: 
def properties(Polygon):
    ''' 
    Computes area, perimeter and diameter of a polygon
    
    Parameter
    ---------
    
    Polygon:  shapely polygon
        A closed polygon 
        
    Returns
    -------
    
    list
        contains area, perimeter and diameter
    '''
    A = Polygon.area
    Perimeter = Polygon.length
    Diameter = maximal_diameter(Polygon)
    return [A, Perimeter, Diameter]

    
def classification(OldFileName, NewFileName):
    """
    Classification of polygons in OldFileName which represent the gaps which 
    should be filled. All polygons in one class get the same kriging radius.
    The radius is saved as an attribute in the output file NewFileName.
    
    Parameter
    ---------
    OldFileName  :  str
        Path of the input .shp file
    NewFileName  : str
        Path of the output .shp file
    """
    with fiona.open(OldFileName, 'r') as Gaps:
        NumPolygons = len(list(Gaps))
        Properties = np.empty([NumPolygons, 3]) # Characteristics of gaps
        Rad = np.empty([NumPolygons, 1]) # Kriging radius
        i = 0
        for Polygon in Gaps:
            P = shape(Polygon['geometry'])
            Properties[i] = properties(P)
            Rad[i] = maximal_diameter(P)+500
            i = i + 1
    Properties = np.divide(Properties, np.absolute(Properties).max(0)) 
        # scaling property values
    NumClusters = 40
    if NumClusters > NumPolygons:
        NumClusters = NumPolygons
    Cluster = KMeans(n_clusters=NumClusters)
    Cluster.fit(Properties)
    Classes = Cluster.predict(Properties)
    for j in range(NumClusters):
        if Rad[Classes == j].shape[0]>1:
            Rad[Classes == j] = np.max(Rad[Classes == j])
    i = 0
    with fiona.open(OldFileName, 'r') as Gaps:
        Schema = Gaps.schema
        Schema['properties'].update({u'KrigRad': 'int'})
        with fiona.open(NewFileName, 'w', 'ESRI Shapefile', Schema, Gaps.crs)\
                as NewGaps:
            for Polygon in Gaps:
                Polygon['properties'].update({u'KrigRad': int(Rad[i])})
                i = i + 1
                NewGaps.write({'properties': Polygon['properties'],
                           'geometry': Polygon['geometry']})


def extent(FileName):
    """
    Computation of the extent of a raster file (GeoTiff).
    
    Parameter
    ---------
    
    FileName  : str
        Path of the inout file
        
    Returns
    -------
    
    list
        list of the coordinates of the corners of the GeoTiff
    """
    Data = gdal.Open(FileName, gdal.GA_ReadOnly)
    GeoT = Data.GetGeoTransform()
    MinX = GeoT[0]
    MaxY = GeoT[3]
    MaxX = MinX + GeoT[1] * Data.RasterXSize
    MinY = MaxY + GeoT[5] * Data.RasterYSize
    Data = None
    return [MinX, MinY, MaxX, MaxY]
    
 
def main():
    print_time("Computation started")
    Platform = platform.platform()
    # unknown platform results in 
    #   ERROR 4: : No such file or directory
    #   Traceback (most recent call last):
    #   File "int.py", line 392, in <module>
    #   main()
    #   File "int.py", line 359, in main
    #   PriorFile = choose_prior(InFile, InputEPSG, Platform)
    #   File "int.py", line 93, in choose_prior
    #   SRS.ImportFromWkt(DataSet.GetProjectionRef())
    #   AttributeError: 'NoneType' object has no attribute 'GetProjectionRef'
    if Platform == "Windows-10-10.0.15063-SP0":
        Platform = "Windows"
    elif Platform == "Darwin-17.2.0-x86_64-i386-64bit":
        Platform = "Mac"
    elif "Linux" in Platform:
        Platform = "Ollie"
    else:
        Platform = "unknown"
    Resolution = 250
    InFile = sys.argv[1]
    OutFile = sys.argv[2]
    DirName = OutFile.rpartition("/")[0]
    BaseName = InFile.rpartition("/")[2]
    BaseName = BaseName.rpartition(".")[0]
    PriorFile=DirName+"/rignot_vx_rand.tif"
    os.system("gdalwarp -q -overwrite -dstnodata nan -tr "+str(Resolution)+" "+str(Resolution)+" "+DirName+"/prior.tif")
    PriorFile=DirName+"/prior.tif"
    Extent = extent(PriorFile)
    Extent = str(Extent[0])+" "+str(Extent[1])+" "+str(Extent[2])+" "+str(Extent[3])
    # TODO: remove EPSG code # -t_srs EPSG:3413
    os.system("gdalwarp -overwrite -t_srs EPSG:3031 -dstnodata nan -tr "+str(Resolution)+" "+str(Resolution)+" -te "+Extent+" "+InFile+" "+DirName+"/"+BaseName+"_input.tif")
    InFile = DirName+"/"+BaseName+"_input.tif"
    # TODO: InputProj4 not needed
    Input, InputGeoT, InputWKT, InputEPSG, InputProj4 = read(InFile, GetGeoT=True, GetWKT=True, GetEPSG=True, GetProj4=True)
    #PriorFile = choose_prior(InFile, InputEPSG, Platform)
    #os.system("gdalwarp -q -overwrite -dstnodata nan -tr "+str(Resolution)+" "+str(Resolution)+" -te "+Extent+" "+PriorFile+" "+DirName+"/prior.tif")
    PriorFile=DirName+"/prior.tif"
    Prior = read(PriorFile)[0]
    Input = combine(Input, Prior)
    #print_time("Add edges")
    #Input = add_edge(Input, Prior)
    #print(Input[0,0])
    #write(Input, InFile, InputGeoT, InputWKT)
    print_time("Vectorization")
    Isnan = np.isnan(Input)
    IsnanFile = DirName+"/isnan.tif"
    write(Isnan, IsnanFile, InputGeoT, InputWKT)
    VecFile = DirName+"/vectorized"
    if not os.path.exists(VecFile):
        os.makedirs(VecFile)
    os.system("./gdal_polygonize.py -q "+IsnanFile+" "+VecFile+" -f \"ESRI Shapefile\" ")
    print_time("Computation of kriging radii")
    RadiusSHPFile = DirName+"/kriging_radius"
    RadiusFile = DirName+"/kriging_radius.tif"
    classification(VecFile, RadiusSHPFile)
    os.system("gdal_rasterize -q -a KrigRad -l kriging_radius -a_srs EPSG:"+InputEPSG+" -tr "+str(Resolution)+" "+str(Resolution)+" -te "+Extent+" "+RadiusSHPFile+"/kriging_radius.shp "+RadiusFile)
    print_time("Start Kriging")
    #TODO: Add kriging command for other systems
    if Platform == "Windows":
        Command = ["C:/PROGRA~1/R/R-34~1.2/bin/x64/Rscript.exe", "kriging_windows.R", InFile, RadiusFile, OutFile, DirName]
        Out=subprocess.check_output(Command)
        print(Out.decode())
    elif Platform == "Mac":
        os.system("./kriging_mac.R "+InFile+" "+RadiusFile+" "+OutFile+" "+DirName)
    elif Platform == "Ollie":
        os.system("Rscript kriging.R "+DirName+"/"+BaseName+"_input.tif"+" "+DirName+"/kriging_radius.tif"+" "+DirName+"/output.tif"+" "+DirName)
    else:
        print("Not ready for this system")
        #os.system("./kriging.R "+InFile+" "+RadiusFile+" "+OutFile+" "+DirName)

    
if __name__ == "__main__":
    main()
