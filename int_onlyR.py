import os
import sys

InFile = sys.argv[1]
OutFile = sys.argv[2]
DirName = OutFile.rpartition("/")[0]
BaseName = InFile.rpartition("/")[2]
BaseName = BaseName.rpartition(".")[0]
InFile = DirName+"/"+BaseName+"_input.tif"
RadiusFile = DirName+"/kriging_radius.tif"

#os.system("gdal_translate -of XYZ "+InFile+" "+DirName+"/"+BaseName+"_input.xyz")
#os.system("gdal_translate -of XYZ "+RadiusFile+" "+DirName+"/kriging_radius.xyz")
os.system("Rscript kriging.R "+DirName+"/"+BaseName+"_input.tif"+" "+DirName+"/kriging_radius.tif"+" "+DirName+"/output.tif"+" "+DirName)
#os.system("gdal_translate "+DirName+"/output.tif "+OutFile)
