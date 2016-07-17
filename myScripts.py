#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################
# AUTHOR: Mateusz Kędzior
# PURPOSE: Python scripts to perform processing of Sentinel-1 data using ESA SNAP
# PREREQUISITES: 
# - install ESA SNAP, go to terminal and type:
#    cd ~/.snap/snap-python/snappy
#    sudo /usr/bin/python setup.py install
# - install HDF5 libraries for Java:
#    sudo apt install libjhdf5-jni libjhdf5-java
# - java_max_mem setting in ~/.snap/snap-python/snappy/snappy.ini is not interpreted by snappy
#   so I set _JAVA_OPTIONS in the first lines of scripts to use 4 GB of RAM

# DO NOT forget that snappy for ESA SNAP is not Google library with exactly the same name!!
# Dokumentacja API do SNAP:
# http://step.esa.int/docs/v3.0/apidoc/desktop/
#############################################################

import sys, os

# To avoid RuntimeError: java.lang.OutOfMemoryError: Java heap space
print("Current _JAVA_OPTIONS: '" + os.environ.get('_JAVA_OPTIONS','Not Set'))
print("will be changed to '-Xmx4096m' to avoid OutOfMemoryError")
os.environ["_JAVA_OPTIONS"] = "-Xmx4096m"
os.system('export _JAVA_OPTIONS=-Xmx4096m')

import snappy
from snappy import ProductIO
from snappy import GPF
from snappy import jpy

from os.path import expanduser
home = expanduser("~")

SentinelPath = os.path.join(home,"Testy")
SentinelFile = os.path.join(SentinelPath, "S1A_IW_GRDH_1SDV_20160512T161044_20160512T161109_011228_010FA8_C584.zip")

def newFilepath(Filepath, prefix):
	return os.path.join(os.path.dirname(Filepath), "_".join([prefix,os.path.basename(Filepath)[0:45]]) + ".beam")

def getDateFromSMOSfileName(SMOSfile1):
	import re, os
	result = re.findall("CLF3.._(\d{8}).*",os.path.basename(SMOSfile1))
	if (len(result) !=1):
		print "Unable to get date from SMOS file name: " + os.path.basename(SMOSfile1)
		return False
	else:
		return result[0]

def getNewSMOSfileName(SMOSfile1, SMOSfile2,destination):
	import os
	date1 = getDateFromSMOSfileName(SMOSfile1)
	date2 = getDateFromSMOSfileName(SMOSfile2)
	return os.path.join(destination, "_".join(["SMOS",date1,date2]) + ".beam")

def getSigma(SentinelFile):
	# calculate sigma (radar backscatter)
	# in ESA SNAP desktop: Radar --> Radiometric --> Calibrate
	if os.path.exists(SentinelFile):
		# Read sourceProduct
		sourceProduct = snappy.ProductIO.readProduct(SentinelFile)
		# Use calibration operator - I've taken "org.esa.s1tbx.calibration.gpf.CalibrationOp" from the help window
		CalibrationOp = jpy.get_type("org.esa.s1tbx.calibration.gpf.CalibrationOp")
		CalOp = CalibrationOp()
		CalOp.setParameterDefaultValues()
		CalOp.setSourceProduct(sourceProduct)
		CalOp.setParameter('doExecute', True)
		# Don't need to create the target product. It is created by the operator.
		targetProduct = CalOp.getTargetProduct()
		newFile=newFilepath(SentinelFile, "calibrated")
		print("Starting writing to the file: " + newFile)
		snappy.ProductIO.writeProduct(targetProduct, newFile, 'BEAM-DIMAP')
		return newFile

def getSubset(SentinelFilePath):
	###################################################
	## Configuration
	#Polygon should describe part of the Eastern Poland
	wkt = "POLYGON((23.00 52.00,24.00 52.00,24.00 52.25,23.00 52.25,23.00 52))"
	#Prefix added to new file:
	prefix = "SUBSET"
	# type of output:
	OutputType = [".dim","BEAM-DIMAP"]
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
	newFile=newFilepath(SentinelFile, prefix)
	if os.path.exists(newFile + OutputType[0]):
		print("It seems that subset of your data already exists. Bye!")
		return
	else:
		print("Starting writing to the file: " + newFile)
		ProductIO.writeProduct(sub_product, newFile, OutputType[1])
		return newFile + OutputType[0]

def getDiff(file1,file2,destination, band='Soil_Moisture'):
	import snappy
	from snappy import GPF
	from snappy import ProductIO

	products = [ snappy.ProductIO.readProduct(file1), snappy.ProductIO.readProduct(file2) ]
	#verify if products contain selected band
	for prod in products:
		if not (band in prod.getBandNames()):
			print(band + " not in the " + prod.getName() + " Exiting")
			return

	GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

	HashMap = jpy.get_type('java.util.HashMap')
	BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

	targetBand1 = BandDescriptor()
	targetBand1.name = 'Difference'
	targetBand1.type = 'float32'
	##  index is zero-based, so index 1 refers to the second product
	targetBand1.expression = band + ' - $sourceProduct.1.' + band

	targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
	targetBands[0] = targetBand1

	parameters = HashMap()
	parameters.put('targetBands', targetBands)

	## More at http://forum.step.esa.int/t/calculate-the-difference-or-division-between-bands-in-two-different-products
	result = GPF.createProduct('BandMaths', parameters, products)

	resultFile = getNewSMOSfileName(file1, file2,destination)
	ProductIO.writeProduct(result, resultFile, 'BEAM-DIMAP')
	return resultFile

def getDivision(file1,file2,band='Soil_Moisture'):
    pass
