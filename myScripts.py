#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################
# AUTHOR: Mateusz Kędzior
# PURPOSE: Python scripts to perform Sentinel-1 processing data using ESA SNAP
# PREREQUISITES:
# - install ESA SNAP, go to terminal and type:
#    cd ~/.snap/snap-python/snappy
#    sudo /usr/bin/python setup.py install
# - install HDF5 libraries for Java:
#    sudo apt install libjhdf5-jni libjhdf5-java
# - java_max_mem setting in ~/.snap/snap-python/snappy/snappy.ini
#   is not interpreted by snappy
#   so I set _JAVA_OPTIONS in the first lines of scripts to use 4 GB of RAM
# - to avoid errors jhdf5 errors as described here: http://forum.step.esa.int/t/snappy-hdf5-error/867/3
#    execute following lines:
#    SNAP_HOME=~/snap
#    cd $SNAP_HOME/snap/modules/lib/x86_64
#    ln -s ../amd64/libjhdf.so
#    ln -s ../amd64/libjhdf5.so

# DO NOT forget that snappy for ESA SNAP is not Google library!!
# Dokumentacja API do SNAP:
# http://step.esa.int/docs/v3.0/apidoc/desktop/
#############################################################

import os
import sys
reload(sys)  
sys.setdefaultencoding('utf8')

# Set below-normal priority, so that computer remain responsive during computations
# You can check how to do that on non-Unix like machines at:
# http://stackoverflow.com/questions/1023038/change-process-priority-in-python-cross-platform
os.nice(10)

# To avoid RuntimeError: java.lang.OutOfMemoryError: Java heap space
print(("Current _JAVA_OPTIONS: '" + os.environ.get('_JAVA_OPTIONS', 'Not Set')))
print("will be changed to '-Xmx4096m' to avoid OutOfMemoryError")
os.environ["_JAVA_OPTIONS"] = "-Xmx4096m"
os.system('export _JAVA_OPTIONS=-Xmx4096m')

##### Global parameters:
# output files format:
OutputType = [".dim", "BEAM-DIMAP"]
# Area - Polygon should describe part of the Eastern Poland
wkt = "POLYGON((23.00 52.00,24.00 52.00,24.00 52.25,23.00 52.25,23.00 52))"
# prefixes added to file names:
prefixes = ["calibrated", "subset"]
# pixel spacing
destinationPS = float(500)
## SMOS pixel size: Pixel Size = (0.260303676128387,-0.314965192009421)
SMOSPS = 28963
SentinelPS = 10.0

import snappy
from snappy import ProductIO
#from snappy import GPF
from snappy import jpy

# Sample file used in testing:
from os.path import expanduser
home = expanduser("~")
SentinelPath = os.path.join(home, "Testy")
SentinelFile = os.path.join(SentinelPath,
"S1A_IW_GRDH_1SDV_20160512T161044_20160512T161109_011228_010FA8_C584.zip")
smallFile = os.path.join(SentinelPath,"collocation/CLF33A_20160514_collocation_20160517_.data/Soil_Moisture_M.img")
sampleDimFile = os.path.join(SentinelPath,"CLF33A_20160514_collocation_20160517_.dim")


def getAllowedFormats():
	# Allowed formats to write: GeoTIFF-BigTIFF,HDF5,Snaphu,BEAM-DIMAP,
	# GeoTIFF+XML,PolSARPro,NetCDF-CF,NetCDF-BEAM,ENVI,JP2,
	# Generic Binary BSQ,Gamma,CSV,NetCDF4-CF,GeoTIFF,NetCDF4-BEAM
	ProductIOPlugInManager = jpy.get_type(
		'org.esa.snap.core.dataio.ProductIOPlugInManager')

	read_plugins = ProductIOPlugInManager.getInstance().getAllReaderPlugIns()
	write_plugins = ProductIOPlugInManager.getInstance().getAllWriterPlugIns()

	writerFormats = ""
	readerFormats = ""

	while write_plugins.hasNext():
		plugin = next(write_plugins)
		writerFormats = writerFormats + plugin.getFormatNames()[0] + ","

	while read_plugins.hasNext():
		plugin = next(read_plugins)
		readerFormats = readerFormats + plugin.getFormatNames()[0] + ","

	print(("Allowed formats to write: " + writerFormats))
	print(("Allowed formats to read: " + readerFormats))


def newFilepath(Filepath, prefix, limited=True):
	directory = os.path.join(os.path.dirname(Filepath),prefix)
	if not os.path.exists(directory):
		os.makedirs(directory)
	baseName = os.path.basename(Filepath)[0:45] if limited else os.path.splitext(os.path.basename(Filepath))[0]
	return os.path.join(directory,
	"_".join([prefix, baseName]) + OutputType[0])

def getDateFromFileName(FileName):
	import re
	daty = re.findall('(\d{8})',FileName)
	if daty == []:
		daty = re.findall('(\d{4}.\d{2}.\d{2})',FileName)
		return daty[0].replace(".","")
	else:
		if validateDate(daty[0]):
			return daty[0]
		else:
			daty = re.findall('(\d{4}.\d{2}.\d{2})',FileName)
			return daty[0].replace(".","")

def validateDate(date_text, dateFormat='%Y%m%d'):
	import datetime
	try:
		datetime.datetime.strptime(date_text, dateFormat)
		return True
	except ValueError:
		return False #raise ValueError("Incorrect data format, should be: "+ dateFormat)

def getListOfFiles(folderPath):
	import os
	return os.listdir(folderPath)

def getListOfDatesFromFileName(folderPath):
	mm = getListOfFiles(folderPath)
	posort = []
	for i in mm:
		posort.append([getDateFromFileName(i),i])
	posort.sort()
	return posort

def getFilterredListOfDatesAndFiles(folderPath,extension=".nc"):
	import os
	myList = getListOfDatesFromFileName(folderPath)
	for j in myList:
		if (os.path.splitext(j[1])[1] == extension):
			print "\t".join([j[0],j[1],os.path.splitext(j[1])[1]])

def writeFilterredListToFile(folderPath,extension=".nc"):
	import os
	charset = ""
	myList = getListOfDatesFromFileName(folderPath)
	f = open(os.path.join(folderPath,"list.txt"),'w')
	for j in myList:
		if (os.path.splitext(j[1])[1] == extension):
			charset = charset + "\t".join([j[0],j[1],os.path.splitext(j[1])[1]]) + "\n"
	f.write(charset)
	f.close()

def unpackAllAndRemoveAllArchives(folderPath,extension="tgz"):
	import os, glob, tarfile
	path = os.path.join(folderPath,"*." + extension)
	FileList = glob.glob(path)
	for archive in FileList:
		tar=tarfile.open(archive)
		tar.extractall(path=folderPath)
		tar.close()
		os.remove(archive)


def getDateFromSMOSfileName(SMOSfile1):
	import re
	import os
	result = re.findall("(20\d{6}).*", os.path.basename(SMOSfile1))
	if (len(result) != 1):
		print(("Unable to get date from SMOS file name: "
		+ os.path.basename(SMOSfile1)))
		return False
	else:
		return result[0]


def getNewFileName(SMOSfile1, SMOSfile2, destination, operation, band, filetype, getFullName=False):
	import os
	if getFullName:
		date1 = os.path.splitext(os.path.basename(SMOSfile1))[0]
		date2 = os.path.splitext(os.path.basename(SMOSfile2))[0]
	else:
		date1 = getDateFromSMOSfileName(SMOSfile1)
		date2 = getDateFromSMOSfileName(SMOSfile2)
	directory = os.path.join(destination,operation)
	if not os.path.exists(directory):
		os.makedirs(directory)
	return os.path.join(directory,
	"_".join([filetype, date1, operation, date2,band]) + OutputType[0])


def writeToLog(message, messageType='ERROR'):
	import logging,os
	log_path = os.path.join(os.path.expanduser("~"),'Dropbox/DyzagregacjaSMOS/logi.log')
	logger = logging.getLogger('Dyzagregacja')
	hdlr = logging.FileHandler(log_path)
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.setLevel(logging.INFO)
	print("Writing to log file: '{0}'".format(log_path))
	if (messageType == 'ERROR'):
		logger.error(message)
	elif (messageType == 'WARNING'):
		logger.warning(message)
	else:
		logger.info(message)
	logger.removeHandler(hdlr)


def convert_bytes(num):
	"""
	this function will convert bytes to MB.... GB... etc
	"""
	for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
		if num < 1024.0:
			return "%3.1f %s" % (num, x)
		num /= 1024.0

def get_size(start_path = '.'):
	# returns total size of folder in bytes
	total_size = 0
	for dirpath, dirnames, filenames in os.walk(start_path):
		for f in filenames:
			fp = os.path.join(dirpath, f)
			total_size += os.path.getsize(fp)
	return total_size

def get_whole_Product_size(file_path):
	total_size = file_size_in_bytes(file_path)
	if (os.path.splitext(file_path)[1] == '.dim'):
		total_size += get_size(get_data_path(file_path))
	if (os.path.splitext(file_path)[1] == '.nc'):
		sFile =(os.path.splitext(file_path)[0]).replace("_1.DBL","").replace("ext-","") + '.HDR'
		total_size += get_size(sFile)
	return convert_bytes(total_size)

def get_data_path(file_path):
	data_path = os.path.splitext(file_path)[0] + '.data'
	if (not os.path.isdir(data_path)):
		message = "There is NO following foler '{0}'. Please ensure where data for '{1}' file are".format(data_path, file_path)
		writeToLog(message)
	return data_path

def file_size_in_bytes(file_path):
	"""
	this function will return the file size
	"""
	if os.path.isfile(file_path):
		file_info = os.stat(file_path)
		return float(file_info.st_size)
	else:
		return float(0)

def file_size(file_path):
	return convert_bytes(file_size_in_bytes(file_path))

def removeProduct(file_path):
	import shutil
	if (os.path.exists(file_path)):
		message = "Trying to remove file '{0}' - file size {1}, whole product size \t {2}".format(os.path.basename(file_path), file_size(file_path), get_whole_Product_size(file_path))
		writeToLog(message,"info")
		if (os.path.splitext(file_path)[1] == '.dim'):
			dirToRem = get_data_path(file_path)
			shutil.rmtree(dirToRem) # will delete a directory and all its contents.
		os.remove(file_path)
	else:
		writeToLog("Trying to remove non-existing file: {0}".format(file_path),"warning")

def createMap(raster, vmax, vmin, output, shapefile=None, title=None):
	###################################################################
	# Author: Mateusz Kędzior
	# Creates image from raster and shapefile
	# Based on: https://gist.github.com/jdherman/7434f7431d1cc0350dbe
	######
	# TODO: Consider rewriting to pyQGIS
	# (http://docs.qgis.org/testing/en/docs/pyqgis_developer_cookbook/composer.html)
	#####
	# Prerequisities:
	# sudo apt-get install python-mpltoolkits.basemap
	##################################################################
	## Sample files for testing (comment in gedit: CTRL + M, uncomment: CTRL + SHIFT + M)
	#import os
	#from os.path import expanduser
	#home = expanduser("~")

	#SMOSfile = os.path.join(home,"Dropbox/Dane SMOS CATDS dla Wisły/DA_TC_MIR_CL_33/EXT-SM_RE02_MIR_CLF33A_20101231T000000_20120102T235959_272_001_7/ext-SM_RE02_MIR_CLF33A_20101231T000000_20120102T235959_272_001_7_1.DBL.nc")
	#SMOSraster = 'NETCDF:"' + SMOSfile + '":Soil_Moisture'
	#SentinelRaster = os.path.join(home,"Testy/calibrated_S1A_IW_GRDH_1SDV_20160512T161044_20160512T161.data/Sigma0_VH.img")
	#SMAPfile = os.path.join(home,"SMOSSMAPAquarius/SMAP/2015.04.15/SMAP_L3_SM_AP_20150415_R13080_001.h5")
	#SMAPraster = 'HDF5:"' + SMAPfile + '"://Soil_Moisture_Retrieval_Data/soil_moisture'

	#Aquariusfile = os.path.join(home,"SMOSSMAPAquarius/Aquarius/2015.04.13/Q2015103.L3m_DAY_SOILM_V4.0_rad_sm_1deg")
	#Aquariusraster = 'HDF5:"' + Aquariusfile + '"://l3m_data'
	#vmin = 0
	#vmax = 1
	#output = os.path.join(home,"testy.png")
	#shapefile = os.path.join(home,"Dropbox/mapy/dorzecze_Wisły")
	#createMap(SMOSraster, vmax, vmin, output, shapefile)
	#createMap(SentinelRaster, vmax, vmin, output)
	###################################################################
	
	# Set font which contains polish characters:
	import matplotlib
	matplotlib.rc('font', family='Arial')

	from osgeo import gdal, osr
	import matplotlib.pyplot as plt
	import numpy as np
	from mpl_toolkits.basemap import Basemap

	# Clear plot
	plt.clf()
	m = None
	cmap = None
	im = None

	
	# By default, osgeo.gdal returns None on error, and does not normally raise informative exceptions
	gdal.UseExceptions()
	
	gdata = gdal.Open(raster)
	geo = gdata.GetGeoTransform()

	factor = float(gdata.GetMetadataItem('Soil_Moisture#scale_factor')) if gdata.GetMetadataItem('Soil_Moisture#scale_factor') != None else float(1)

	xres = geo[1]
	yres = geo[5]

	# Get "natural" block size, and total raster XY size. 
	band = gdata.GetRasterBand(1)
	block_sizes = band.GetBlockSize()
	x_block_size = block_sizes[0]
	y_block_size = block_sizes[1]
	xsize = band.XSize
	ysize = band.YSize
	print('x_block_size: {0}, y_block_size: {1}.'.format(x_block_size, y_block_size))
	print('xsize: {0}, ysize: {1}.'.format(xsize, ysize))

	if (xsize < 5000):
		data = gdata.ReadAsArray()
		data = data * factor
		print('Whole data min: {0}, max: {1}, mean: {2}.'.format(data.min(), data.max(), data.mean()))
	else:
		#########################################################
		## TODO: for big rasters such as Sentinel-1:
		## Solution adapted from http://gis.stackexchange.com/questions/211611/python-gdal-handling-big-rasters-and-avoid-memoryerrors
		## It seems that I still receive Memory Error 
		y_block_size_NEW = int(round(y_block_size/200)) if y_block_size > 200 else y_block_size
		x_block_size_NEW = int(round(x_block_size/200)) if x_block_size > 200 else x_block_size

		# Create temporal raster
		raster_srs = osr.SpatialReference()
		raster_srs.ImportFromWkt(gdata.GetProjectionRef())

		format = "GTiff"
		driver = gdal.GetDriverByName( format )
		# TODO: seems that I should create smaller temporal raster (?)
		dst_ds = driver.Create("original_blocks.tif", xsize, ysize, 1, band.DataType )

		dst_ds.SetGeoTransform(geo)
		dst_ds.SetProjection(raster_srs.ExportToWkt())

		blocks = 0 
		for y in xrange(0, ysize, y_block_size_NEW):
			#print blocks
			if y + y_block_size_NEW < ysize:
				rows = y_block_size_NEW
			else:
				rows = ysize - y
			for x in xrange(0, xsize, x_block_size_NEW):
				if x + x_block_size_NEW < xsize:
					cols = x_block_size_NEW	
				else:
					cols = xsize - x
				# Seems that some kind of average should be calculated here
				array = band.ReadAsArray(x, y, cols, rows)
				try:
					array[array>0]=1
					#print "we got them"
				except:
					print "could not find them"
				dst_ds.GetRasterBand(1).WriteArray(array, x, y)
				del array
				blocks += 1

		data = dst_ds.ReadAsArray() * factor
		# TODO: Remove temporal raster?
		#########################################################
	#23.00 52.00,24.00 52.00,24.00 52.25,23.00 52.25,23.00 52
	m = Basemap(llcrnrlon=22.75,llcrnrlat=51.75,urcrnrlon=24.25,urcrnrlat=52.25)
	#m = Basemap(llcrnrlon=17.00,llcrnrlat=48.75,urcrnrlon=25.25,urcrnrlat=54.50)

	if shapefile is not None:
		m.readshapefile(shapefile,'shp',drawbounds=True, color='0.3')
	xmin = geo[0] + xres * 0.5
	xmax = geo[0] + (xres * gdata.RasterXSize) - xres * 0.5
	ymin = geo[3] + (yres * gdata.RasterYSize) + yres * 0.5
	ymax = geo[3] - yres * 0.5
	x,y = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
	x,y = m(x,y)
	cmap = plt.cm.gist_rainbow
	cmap.set_under ('1.0')
	cmap.set_bad('0.8')
	im = m.pcolormesh(x,y, data.T, cmap=cmap, vmin=vmin, vmax=vmax)
	cb = plt.colorbar( orientation='vertical', fraction=0.10, shrink=0.7)
	if title is not None:
		plt.title(title)
	plt.savefig(output) # to take less space add: bbox_inches='tight', pad_inches=0

def createMAPsForFolder(path, fileMask, outputPath, fileName, whatADD=[], shapefile=None):
	import fnmatch
	import os

	vmax = 1
	vmin = 0

	matches = []
	for root, dirnames, filenames in os.walk(path):
		for filename in fnmatch.filter(filenames, fileMask):
			matches.append(os.path.join(root, filename))

	for myFile in matches:
		if whatADD != []:
			raster = whatADD[0] + myFile + whatADD[1]
		else:
			raster = myFile
		print('Trying to get date for following file name: {0}'.format(myFile))
		myDate = getDateFromFileName(myFile)
		output = os.path.join(outputPath,myDate + fileName)
		print('raster: {0}, output: {1}.'.format(raster, output))
		createMap(raster, vmax, vmin, output, shapefile)

def getSigma(SentinelFile):
	# calculate sigma (radar backscatter)
	# in ESA SNAP desktop: Radar --> Radiometric --> Calibrate
	if os.path.exists(SentinelFile):
		newFile = newFilepath(SentinelFile, prefixes[0])
		if (not os.path.exists(newFile)):
			# Read sourceProduct
			sourceProduct = snappy.ProductIO.readProduct(SentinelFile)
			# Use calibration operator - I've taken:
			# "org.esa.s1tbx.calibration.gpf.CalibrationOp" from the help window
			CalibrationOp = jpy.get_type("org.esa.s1tbx.calibration.gpf.CalibrationOp")
			CalOp = CalibrationOp()
			CalOp.setParameterDefaultValues()
			CalOp.setSourceProduct(sourceProduct)
			CalOp.setParameter('doExecute', True)
			# Don't need to create the target product. It is created by the operator.
			targetProduct = CalOp.getTargetProduct()
			print(("Starting writing to the file: " + newFile))
			snappy.ProductIO.writeProduct(targetProduct, newFile, OutputType[1])
			sourceProduct.dispose()
			targetProduct.dispose()
		else:
			print(("File already exists. Exit without changes."))
		return newFile


def getSubset(SentinelFile):
	#Initialize:
	print(("Please exceute getSubset method AFTER executing getSigma (after using Calibration)"))
	SubsetOp = snappy.jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
	WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader')
	geom = WKTReader().read(wkt)
	op = SubsetOp()
	# read source product and set properties:
	product = ProductIO.readProduct(SentinelFile)
	op.setSourceProduct(product)
	op.setGeoRegion(geom)
	sub_product = op.getTargetProduct()
	# Ensure that file does not exists:
	newFile = newFilepath(SentinelFile, prefixes[1])
	if os.path.exists(newFile):
		print("It seems that subset of your data already exists. Bye!")
	else:
		print(("Starting writing to the file: " + newFile))
		ProductIO.writeProduct(sub_product, newFile, OutputType[1])
	product.dispose()
	sub_product.dispose()
	return newFile


def getOperation(file1, file2, destination, operation, band=['Soil_Moisture','Soil_Moisture']):
	import snappy
	from snappy import GPF
	from snappy import ProductIO

	products = [snappy.ProductIO.readProduct(file1),
	snappy.ProductIO.readProduct(file2)]
	#verify if products contain selected band
	if not (band[0] in products[0].getBandNames()):
		writeToLog(band[0] + " not in the " + products[0].getName() + " Exiting")
		writeToLog("Available bands: {0} in file {1}".format(products[0].getBandNames(),os.path.basename(file1)))
		return
	if not (band[1] in products[1].getBandNames()):
		writeToLog(band[1] + " not in the " + products[1].getName() + " Exiting")
		writeToLog("Available bands: {0} in file {1}".format(products[1].getBandNames(),os.path.basename(file2)))
		return

	GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

	HashMap = jpy.get_type('java.util.HashMap')
	BandDescriptor = jpy.get_type(
	'org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

	targetBand1 = BandDescriptor()
	targetBand1.name = operation[1]
	targetBand1.type = 'float32'

	if (getProductInfo(file1) == getProductInfo(file2)):
		##  index is zero-based, so index 1 refers to the second product
		expr = "".join([band[0], ' ', operation[0], ' $sourceProduct.1.', band[1]])
	else:
		# in this case we need collocated product
		# first remove old products:
		for prod in products:
			prod.dispose()
		collocated = getCollocated(file1, file2, destination)
		products = [snappy.ProductIO.readProduct(collocated)]
		# Sample expression: 'Sigma0_VH_M - Sigma0_VH_S'
		expr = "{0}_M {1} {2}_S".format(band[0], operation[0], band[1])
	prodlist = ""
	for prod in products:
		prodlist = prodlist + "'{0}'".format(prod.getName())
	writeToLog("getOperation: I will try to calculate using following expression: '{0}' on following products: {1}".format(expr,prodlist), "info")
	targetBand1.expression = expr

	targetBands = jpy.array(
	'org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
	targetBands[0] = targetBand1

	parameters = HashMap()
	parameters.put('targetBands', targetBands)

	## More at http://forum.step.esa.int/t/calculate-the-difference-or-division
	## -between-bands-in-two-different-products
	result = GPF.createProduct('BandMaths', parameters, products)
	
	# TODO: this should be handled in smarter way!!!
	filetype = os.path.basename(file1).split("_")[3]
	resultFile = getNewFileName(file1, file2, destination, operation[1], band[0], filetype)
	ProductIO.writeProduct(result, resultFile, OutputType[1])
	for prod in products:
		prod.dispose()
	result.dispose()
	return resultFile

def getProductInfo(file1):
	import snappy
	from snappy import GPF
	from snappy import ProductIO

	prod = snappy.ProductIO.readProduct(file1)
	bandNames=''
	for i in prod.getBandNames():
		bandNames += "'{0}'".format(i)
	firstBand=prod.getBands()[0]
	width = firstBand.getRasterWidth()
	height = firstBand.getRasterHeight()
	prod.dispose()
	return "Bands: {0}, width = {1}, height = {2}".format(bandNames,width, height)

def getBandFromProduct(file1, bandNumber):
	import snappy
	from snappy import GPF
	from snappy import ProductIO
	# TODO: When I try to read 'Soil_Moisture_M.img' directly (not hdr file), I receive a NoneType object
	prod = snappy.ProductIO.readProduct(file1)
	
	if (len(prod.getBands()) >= bandNumber):
		Band = prod.getBands()[bandNumber]
	else:
		print("Illegal band number")
		Band = None
	# If 'prod.dispose()' line (below) is not commented, I receive an error message when usign this function in getBandStats
	# RuntimeError: java.lang.IllegalArgumentException: The name to be externalized must at least contain one character
	#prod.dispose()
	return Band

def getAllBandsStats(file1, pathToSaveStats=None):
	import snappy
	# TODO: When I try to read 'Soil_Moisture_M.img' directly (not hdr file), I receive a NoneType object
	prod = snappy.ProductIO.readProduct(file1)
	if (not prod):
		errormsg = "getAllBandsStats - Error when reading '{0}' file".format(file1)
		print(errormsg)
		return errormsg
	numberOfBands = len(prod.getBands())
	prodName = prod.getName()
	fileName = prod.getFileLocation().getName()
	if (not pathToSaveStats):
		pathToSaveStats = os.path.join(SentinelPath,"BandStatistics.csv")
	####
	for bandNumber in range(numberOfBands):
		Band = prod.getBands()[bandNumber]
		stats = Band.getStx()
		message = "FileName,Product,BandName,Min,Max,Avg,StdDev,CV,NumberOfPixels,TotalNumberOfPixels:\t" + ("\t".join(["{0}","{1}","{2}","{3}","{4}","{5}","{6}","{7}","{8}", "{9}"])).format(fileName, prodName, Band.getName(), stats.getMinimum(), stats.getMaximum(), stats.getMedian(), stats.getStandardDeviation(), stats.getCoefficientOfVariation(),int(Band.getNumDataElems()), int(stats.getHistogram().getTotals()[0]))
		with open(pathToSaveStats, "a") as myfile:
			myfile.write(message)
			myfile.write("\n")
	###
	print("Stats saved in '{0}'".format(pathToSaveStats))
	prod.dispose()
	return "Statistics for all '{0}' bands of product '{1}' has been saved in '{2}'".format(numberOfBands,prodName,pathToSaveStats)

def getBandHistogram(file1, bandNumber = 0):
	Band = getBandFromProduct(file1, bandNumber)
	stats = Band.getStx()
	return stats.getHistogram()

def get_envi_header_dict(hdr):
	# Function from: http://gis.stackexchange.com/questions/48618/how-to-read-write-envi-metadata-using-gdal
	import re
	#Get all "key = {val}" type matches
	regex=re.compile(r'^(.+?)\s*=\s*({\s*.*?\n*.*?})$',re.M|re.I)
	matches=regex.findall(hdr)

	#Remove them from the header
	subhdr=regex.sub('',hdr)

	#Get all "key = val" type matches
	regex=re.compile(r'^(.+?)\s*=\s*(.*?)$',re.M|re.I)
	matches.extend(regex.findall(subhdr))

	return dict(matches)

def read_envi_hdr(hdr_file):
	if os.path.exists(hdr_file):
		with open(hdr_file, 'r') as content_file:
			content = content_file.read()
	return get_envi_header_dict(content)

def getMetadataValueFromHdr(hdr_file, HDRkey = 'data gain values'):
	metadata = read_envi_hdr(hdr_file)
	if metadata:
		value = metadata.get(HDRkey)
		return value.replace("{","").replace("}","").strip() if value else None
	else:
		return None

def saveHistForFiles(file1, xtitle="Values", ytitle="Probability", title="Band: ", suffix="eng",directorySuffix = None):
	# This just executes 'saveHistogramForFile' function below
	import glob
	if (os.path.splitext(file1)[1] == '.dim'):
		searchIn = os.path.join(get_data_path(file1),"*.img")
		for myFile in glob.glob(searchIn):
			saveHistogramForFile(myFile, xtitle, ytitle, title, suffix,directorySuffix)
	else:
		saveHistogramForFile(file1, xtitle, ytitle, title, suffix,directorySuffix)

def saveHistogramForFile(file1, xtitle="Values", ytitle="Probability", title=None, suffix="eng",directorySuffix = None):
	# LIMITATIONS: This is *not* working with .dim files
	# Sample usage:
	# saveHistogramForFile(smallFile)
	# saveHistogramForFile(smallFile, "Wartości", "Prawdopodobieństwo", "Pasmo: ", "pl")
	from osgeo import gdal
	import numpy as np

	# Set font which contains polish characters:
	import matplotlib
	matplotlib.rc('font', family='Arial')

	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt

	dataset = gdal.Open(file1)
	band = dataset.GetRasterBand(1)
	data = np.squeeze(band.ReadAsArray())
	allElem = np.prod(data.shape)
	NaNElem = np.count_nonzero((np.isnan(data)))
	NaNprcnt = NaNElem/allElem
	data = data[np.logical_not(np.isnan(data))]

	if (os.path.splitext(file1)[1] == '.img'):
		hdrFile = os.path.splitext(file1)[0] + ".hdr"
		value = getMetadataValueFromHdr(hdrFile, 'data gain values')
		data = data * float(value) if value else data
		if (not title == None):
			bandName = getMetadataValueFromHdr(hdrFile, 'band names')
			title = title + bandName if bandName else title

	# TODO: data (scaling) factors should be handled in separate method
	if (file1.startswith("NETCDF")):
		factor = float(dataset.GetMetadataItem('Soil_Moisture#scale_factor')) if dataset.GetMetadataItem('Soil_Moisture#scale_factor') != None else float(1)
		data = data * factor

	# the histogram of the data
	n, bins, patches = plt.hist(data, facecolor='green') #, 10, normed=1, facecolor='green', alpha=0.75)

	if (float(NaNprcnt) > float(0)):
		xtitle = xtitle + ' (NA: {:.2%})'.format(NaNprcnt)
	plt.xlabel(xtitle)
	plt.ylabel(ytitle)
	if (not title == None):
		plt.title(title)
	plt.grid(True)
	NewFileName = os.path.split(os.path.split(file1)[0])[1] + os.path.basename(file1) + "_hist_" + suffix + ".svg"
	directory = os.path.join(SentinelPath,"histograms")
	if (not directorySuffix == None):
		directory = os.path.join(directory,directorySuffix)
	if not os.path.exists(directory):
		os.makedirs(directory)
	NewFullPath = os.path.join(directory,NewFileName)
	plt.savefig(NewFullPath)
	#plt.savefig(os.path.splitext(NewFullPath)[0] + '.pdf')
	plt.clf()


def getCollocated(file1, file2, destination):
	import snappy
	from snappy import GPF
	from snappy import ProductIO

	# TODO: this should be handled in smarter way!!!
	filetype = os.path.basename(file1).split("_")[3]
	destinationPath = getNewFileName(file1, file2, destination, "collocation", "", filetype,True)

	if (not os.path.exists(destinationPath)):
		products = [snappy.ProductIO.readProduct(file1), snappy.ProductIO.readProduct(file2)]

		HashMap = jpy.get_type('java.util.HashMap')
		GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

		parameters = HashMap()
		sourceProducts = HashMap()
		sourceProducts.put("master", products[0])
		sourceProducts.put("slave", products[1])
		parameters.put('renameMasterComponents', True)
		parameters.put('renameSlaveComponents', True)
		parameters.put('masterComponentPattern', "${ORIGINAL_NAME}_M")
		parameters.put('slaveComponentPattern', "${ORIGINAL_NAME}_S")
		parameters.put('resamplingType', "NEAREST_NEIGHBOUR")
		
		result = GPF.createProduct('Collocate', parameters, sourceProducts)
		ProductIO.writeProduct(result,  destinationPath, 'BEAM-DIMAP')
		
		for prod in products:
			prod.dispose()
		writeToLog("Collocated product saved as '{0}' \t {1}".format(os.path.basename(destinationPath), get_whole_Product_size(destinationPath)),"info")
		writeToLog("Input files: {0}, {1}".format(getProductInfo(file1),getProductInfo(file2)),"info")
		writeToLog("Collocated product: {0}".format(getProductInfo(destinationPath)),"info")
	else:
		print("It seems that destination file '{0}' already exists. Bye!".format(os.path.basename(destinationPath)))
	return destinationPath


def getDiff(file1, file2, destination, band=['Soil_Moisture','Soil_Moisture']):
	# TODO: test output from SNAP desktop and from this file
	return getOperation(file1, file2, destination,
		["-", "diff"], band)


def getDivision(file1, file2, destination, band=['Soil_Moisture','Soil_Moisture']):
	return getOperation(file1, file2, destination,
		["/", "div"], band)

def getSum(file1, file2, destination, band=['Soil_Moisture','Soil_Moisture']):
	return getOperation(file1, file2, destination,
		["+", "sum"], band)

# I will use Sentinel and SMOS data
# I will use two functions: getCoarseResProd (to SMOSPS resolution) or getBetterResProd (to destinationPS resolution)

def getCoarseResProd(file1, destination):
	destinationPath = newFilepath(file1, "aggregated")
	return getResampled(file1, destinationPath, resolution=float(SMOSPS))

def getBetterResProd(file1, destination):
	destinationPath = newFilepath(file1, "interpolated")
	return getResampled(file1, destinationPath, resolution=float(destinationPS))

def getResampled(file1, destinationPath, resolution=destinationPS):
	# TODO: this should be tested!!!
	# More info: http://forum.step.esa.int/t/aggregation-and-interpolation-of-sentinel-products-should-i-use-snappy-or-gdal-tools/2522/3
	import snappy
	from snappy import GPF
	from snappy import ProductIO
	if (not os.path.exists(destinationPath)):
		product = snappy.ProductIO.readProduct(file1)
		HashMap = jpy.get_type('java.util.HashMap')
		GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
		parameters = HashMap()
		parameters.put('sourceProduct', product)
		parameters.put('upsampling', "Bilinear")
		parameters.put('downsampling', "Mean")
		# As I checked in SNAP desktop, 'targetResolution' option is sometimes not available
		# and I need to use targetHeight and targetWidth instead
		# RuntimeError: org.esa.snap.core.gpf.OperatorException: Operator 'ResamplingOp': Value for 'Target resolution' must be of type 'Integer'.
		# So I convert it to Integer
		parameters.put('targetResolution', int(resolution))
		result = GPF.createProduct('Resample', parameters, product)
		ProductIO.writeProduct(result,  destinationPath, 'BEAM-DIMAP')

		product.dispose()
	else:
		print("It seems that destination file '{0}' already exists. Bye!".format(os.path.basename(destinationPath)))
	return destinationPath

def getTerrainCorrected(file1, destinationPath, crs='WGS84(DD)'):
	# According to lveci: "GRD products are not terrain corrected. Due to the nature of SAR acquisitions, in order to get an accurate geocoding you need to account for the geometry of the acquisition"
	# "Over scenes where you have a DEM, you should use range Doppler terrain correction"
	# Radar --> Geometric --> Terrain Correction --> Range-Doppler Terrain Correction
	# Save parameters (saveLatLon, saveDEM) means that additional information (elevation, latLon) will be saved in output file which is not needed for further processing
	# once you downloaded a DEM for an area it stays in the aux folder until you delete it manually. So you won't need to download it again when you process data from the same area.
	import snappy
	from snappy import GPF
	from snappy import ProductIO
	if (not os.path.exists(destinationPath)):
		product = snappy.ProductIO.readProduct(file1)
		
		DEM = "SRTM 3Sec"
		HashMap = jpy.get_type('java.util.HashMap')
		GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
		writeToLog("Terrain correction using {0} file will use pixelSpacingInMeter: {1}".format(DEM,destinationPS),"info")
		parameters = HashMap()
		parameters.put('demName', DEM)
		parameters.put('externalDEMApplyEGM', True)
		parameters.put('demResamplingMethod', "BILINEAR_INTERPOLATION")
		parameters.put('imgResamplingMethod', "BILINEAR_INTERPOLATION")
		parameters.put('pixelSpacingInMeter', float(destinationPS))
		parameters.put('mapProjection', crs)
		parameters.put('nodataValueAtSea', True)
		# This is ONLY for saving DEM within output file - downloaded DEM will be NOT removed from .snap\AuxData\DEMs
		parameters.put('saveDEM', False)
		parameters.put('saveLatLon', False)
		parameters.put('saveIncidenceAngleFromEllipsoid', False)
		parameters.put('saveLocalIncidenceAngle', False)
		parameters.put('saveProjectedLocalIncidenceAngle', False)
	
		result = GPF.createProduct('Terrain-Correction', parameters, product)
		ProductIO.writeProduct(result,  destinationPath, 'BEAM-DIMAP')

		product.dispose()
	else:
		print("It seems that destination file '{0}' already exists. Bye!".format(os.path.basename(destinationPath)))
	return destinationPath

def getReprojected(file1, destinationPath, crs='EPSG:4326'):
	import snappy
	from snappy import GPF
	from snappy import ProductIO
	if (not os.path.exists(destinationPath)):
		product = snappy.ProductIO.readProduct(file1)

		HashMap = jpy.get_type('java.util.HashMap')
		GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

		parameters = HashMap()
		parameters.put('crs', crs)
		parameters.put('resampling', "Nearest")
		result = GPF.createProduct('Reproject', parameters, product)
		ProductIO.writeProduct(result,  destinationPath, 'BEAM-DIMAP')
		product.dispose()
	else:
		print("It seems that destination file '{0}' already exists. Bye!".format(os.path.basename(destinationPath)))

	return destinationPath
