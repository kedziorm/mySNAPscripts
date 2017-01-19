#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Below are the sample global parameters for myScripts.py

########################################################
##### Global parameters:
# Area - sample polygon below should describe part of the Eastern Poland
# Derlo station coordinates: 52.16950    23.36930
myExtent = [23.00, 52.00, 24.00, 52.25]
constBorder = 0.25
wkt = "POLYGON(({0:.2f} {1:.2f},{2:.2f} {1:.2f},{2:.2f} {3:.2f},{0:.2f} {3:.2f},{0:.2f} {4}))".format(myExtent[0], myExtent[1], myExtent[2], myExtent[3], int(myExtent[1]))
createMAPparams = [myExtent[0] - constBorder, myExtent[1] - constBorder, myExtent[2] + constBorder, myExtent[3]]

## SMOS pixel size: Pixel Size = (0.260303676128387,-0.314965192009421)
SMOSPS = 28963
SentinelPS = 10.0
histLabels = ["Values", "Frequency", "Wartości", "Liczebność"]
getTerrainCorrected_DEM = "SRTM 3Sec" # Possible also: "ASTER 1sec GDEM", "SRTM 1Sec Grid", "ACE30"
getTerrainCorrected_demResamplingMethod = "BILINEAR_INTERPOLATION" # Possible also: "NEAREST_NEIGHBOUR", "CUBIC_CONVOLUTION", "BICUBIC_INTERPOLATION"
getTerrainCorrected_imgResamplingMethod = "BILINEAR_INTERPOLATION" # Possible also: "NEAREST_NEIGHBOUR", "CUBIC_CONVOLUTION", "BICUBIC_INTERPOLATION"
getCollocatedResamplingType = "NEAREST_NEIGHBOUR" # Possible also: "BILINEAR_INTERPOLATION", "CUBIC_CONVOLUTION", "BISINC_CONVOLUTION"
getReprojectedResampling = "Nearest"

SentinelPath = os.path.join(home, "Testy")
histogramDirectory = os.path.join(home,"Dropbox/DyzagregacjaSMOS/histograms")
maps = os.path.join(home,"Dropbox/DyzagregacjaSMOS/maps")
########################################################
