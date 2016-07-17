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

# Global parameter type of output:
OutputType = [".dim", "BEAM-DIMAP"]

import snappy
from snappy import ProductIO
#from snappy import GPF
from snappy import jpy

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


def getSigma(SentinelFile):
	# calculate sigma (radar backscatter)
	# in ESA SNAP desktop: Radar --> Radiometric --> Calibrate
	if os.path.exists(SentinelFile):
		newFile = newFilepath(SentinelFile, "calibrated")
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
			newFile = newFilepath(SentinelFile, "calibrated")
			print(("Starting writing to the file: " + newFile))
			snappy.ProductIO.writeProduct(targetProduct, newFile, OutputType[1])
		else:
			print(("File already exists. Exit without changes."))
		return newFile


def getSubset(SentinelFilePath):
	###################################################
	## Configuration
	#Polygon should describe part of the Eastern Poland
	wkt = "POLYGON((23.00 52.00,24.00 52.00,24.00 52.25,23.00 52.25,23.00 52))"
	#Prefix added to new file:
	prefix = "SUBSET"
	#############################################
	#Initialize:
	SubsetOp = snappy.jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
	WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader')
	geom = WKTReader().read(wkt)
	op = SubsetOp()
	# read source product and set properties:
	product = ProductIO.readProduct(SentinelFilePath)
	op.setSourceProduct(product)
	op.setGeoRegion(geom)
	sub_product = op.getTargetProduct()
	# Ensure that file does not exists:
	newFile = newFilepath(SentinelFile, prefix)
	if os.path.exists(newFile + OutputType[0]):
		print("It seems that subset of your data already exists. Bye!")
		return
	else:
		print(("Starting writing to the file: " + newFile))
		ProductIO.writeProduct(sub_product, newFile, OutputType[1])
		return newFile + OutputType[0]


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
	return resultFile


def getDiff(file1, file2, destination, band='Soil_Moisture'):
	# TODO: test output from SNAP desktop and from this file
	return getOperation(file1, file2, destination,
		["-", "diff"], band='Soil_Moisture')


def getDivision(file1, file2, destination, band='Soil_Moisture'):
	return getOperation(file1, file2, destination,
		["/", "div"], band='Soil_Moisture')
