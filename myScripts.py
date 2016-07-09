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
# - modify java_max_mem in ~/.snap/snap-python/snappy/snappy.ini

# DO NOT forget that snappy for ESA SNAP is not Google library with exactly the same name!!
# Dokumentacja API do SNAP:
# http://step.esa.int/docs/v3.0/apidoc/desktop/
#############################################################

import snappy
import sys, os
from snappy import ProductIO
from snappy import GPF
from snappy import jpy

from os.path import expanduser
home = expanduser("~")

SentinelPath = os.path.join(home,"Testy")
SentinelFile = os.path.join(SentinelPath, "S1A_IW_GRDH_1SDV_20160512T161044_20160512T161109_011228_010FA8_C584.zip")

def newFilepath(Filepath, prefix):
	return os.path.join(os.path.dirname(Filepath), "_".join([prefix,os.path.basename(Filepath)[0:45]]) + ".beam")

def getSigma(SentinelFile):
	# calculate sigma (radar backscatter)
	# in ESA SNAP desktop: Radar --> Radiometric --> Calibrate
	if os.path.exists(SentinelFile):
		# Read sourceProduct
		sourceProduct = snappy.ProductIO.readProduct(SentinelFile)
		# Use calibration operator - I've taken "org.esa.s1tbx.calibration.gpf.CalibrationOp" from the help window
		CalibrationOp = jpy.get_type("org.esa.s1tbx.calibration.gpf.CalibrationOp")
		CalOp = CalibrationOp()
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
