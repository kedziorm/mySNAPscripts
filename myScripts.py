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
destinationPS = 1000
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


def newFilepath(Filepath, prefix):
	return os.path.join(os.path.dirname(Filepath),
	"_".join([prefix, os.path.basename(Filepath)[0:45]]) + OutputType[0])

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
	result = re.findall("CLF3.._(\d{8}).*", os.path.basename(SMOSfile1))
	if (len(result) != 1):
		print(("Unable to get date from SMOS file name: "
		+ os.path.basename(SMOSfile1)))
		return False
	else:
		return result[0]


def getNewSMOSfileName(SMOSfile1, SMOSfile2, destination, operation):
	import os
	date1 = getDateFromSMOSfileName(SMOSfile1)
	date2 = getDateFromSMOSfileName(SMOSfile2)
	return os.path.join(destination,
	"_".join(["SMOS", date1, operation, date2]) + OutputType[0])


def writeToLog(message, messageType='ERROR'):
	import logging,os
	logger = logging.getLogger('Dyzagregacja')
	hdlr = logging.FileHandler(os.path.join(os.path.expanduser("~"),'Dropbox/DyzagregacjaSMOS/logi.log'))
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.setLevel(logging.INFO)
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


def file_size(file_path):
	"""
	this function will return the file size
	"""
	if os.path.isfile(file_path):
		file_info = os.stat(file_path)
		return convert_bytes(file_info.st_size)

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

	m = Basemap(llcrnrlon=17.00,llcrnrlat=48.75,urcrnrlon=25.25,urcrnrlat=54.50)

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


def getOperation(file1, file2, destination, operation, band='Soil_Moisture'):
	import snappy
	from snappy import GPF
	from snappy import ProductIO

	products = [snappy.ProductIO.readProduct(file1),
	snappy.ProductIO.readProduct(file2)]
	#verify if products contain selected band
	for prod in products:
		if not (band in prod.getBandNames()):
			print((band + " not in the " + prod.getName() + " Exiting"))
			return

	GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

	HashMap = jpy.get_type('java.util.HashMap')
	BandDescriptor = jpy.get_type(
	'org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

	targetBand1 = BandDescriptor()
	targetBand1.name = operation[1]
	targetBand1.type = 'float32'
	##  index is zero-based, so index 1 refers to the second product
	expr = "".join([band, ' ', operation[0], ' $sourceProduct.1.', band])
	targetBand1.expression = expr

	targetBands = jpy.array(
	'org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
	targetBands[0] = targetBand1

	parameters = HashMap()
	parameters.put('targetBands', targetBands)

	## More at http://forum.step.esa.int/t/calculate-the-difference-or-division
	## -between-bands-in-two-different-products
	result = GPF.createProduct('BandMaths', parameters, products)

	resultFile = getNewSMOSfileName(file1, file2, destination, operation[1])
	ProductIO.writeProduct(result, resultFile, OutputType[1])
	for prod in products:
		prod.dispose()
	result.dispose()
	return resultFile


def getDiff(file1, file2, destination, band='Soil_Moisture'):
	# TODO: test output from SNAP desktop and from this file
	return getOperation(file1, file2, destination,
		["-", "diff"], band)


def getDivision(file1, file2, destination, band='Soil_Moisture'):
	return getOperation(file1, file2, destination,
		["/", "div"], band)

def getSum(file1, file2, destination, band='Soil_Moisture'):
	return getOperation(file1, file2, destination,
		["+", "sum"], band)

# I will use Sentinel and SMOS data
# I will use two functions: getCoarseResProd (to SMOSPS resolution) or getBetterResProd (to destinationPS resolution)

def getCoarseResProd(file1, destinationPath):
	return getResampled(file1, destinationPath, resolution=SMOSPS)

def getBetterResProd(file1, destinationPath):
	return getResampled(file1, destinationPath, resolution=destinationPS)

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
		parameters.put('targetResolution', resolution)
		result = GPF.createProduct('Resample', parameters, product)
		ProductIO.writeProduct(result,  destinationPath, 'BEAM-DIMAP')

		product.dispose()
	else:
		print("It seems that destination file already exists. Bye!")
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

		HashMap = jpy.get_type('java.util.HashMap')
		GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

		parameters = HashMap()
		parameters.put('demName', "SRTM 3Sec")
		parameters.put('externalDEMApplyEGM', True)
		parameters.put('demResamplingMethod', "BILINEAR_INTERPOLATION")
		parameters.put('imgResamplingMethod', "BILINEAR_INTERPOLATION")
		parameters.put('pixelSpacingInMeter', destinationPS)
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
		print("It seems that destination file already exists. Bye!")
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
		print("It seems that destination file already exists. Bye!")

	return destinationPath
