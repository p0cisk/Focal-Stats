"""
/***************************************************************************
 FocalStatsDialog
								 A QGIS plugin
 Raster statistics over a specified neighborhood
							 -------------------
		begin				: 2011-11-22
		copyright			: (C) 2011 by Piotr Pociask
		email				: NA
 ***************************************************************************/

/***************************************************************************
 *																		 *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or	 *
 *   (at your option) any later version.								   *
 *																		 *
 ***************************************************************************/
"""

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *
from qgis.gui import *

import numpy as np
import focal
from osgeo import gdal
from dialogs import *
from time import time


from ui_focalstats import Ui_FocalStats
# create the dialog for zoom to point
class FocalStatsDialog(QDialog):
	def __init__(self):
		QDialog.__init__(self)
		# Set up the user interface from Designer.
		self.ui = Ui_FocalStats()
		self.ui.setupUi(self)
		
		self.ui.cbLayers.addItems(getLayersNames())
		QObject.connect( self.ui.bOK, SIGNAL( "clicked()" ), self.calculateStats )
		QObject.connect( self.ui.bBrowse, SIGNAL( "clicked()" ), self.outFile )
		QObject.connect( self.ui.bOpenRaster, SIGNAL( "clicked()" ), self.openRaster )
		
	def outFile(self):
		"""Open a file save dialog and set the output file path."""
		outFilePath = saveDialog(self)
		if not outFilePath:
			return
		self.ui.eOutput.setText(QString(outFilePath))
		
	def openRaster(self):
		rasterPath = openDialog(self)
		if not rasterPath: return
		self.ui.cbLayers.setEditText(QString(rasterPath))
		
	def createSimpleMask(self):
		sizeX = int(self.ui.cbX.currentText())
		
		if self.ui.useY.isChecked():
			sizeY = int(self.ui.cbY.currentText())
		else:
			sizeY = sizeX
			
		return np.ones((sizeY, sizeX))
		
	def calculateStats(self):
		if self.ui.eOutput.text() == '':
			QMessageBox.information(None, 'Focal statistics', 'Choose output raster file!')
			return
		elif self.ui.cbLayers.currentText() == '':
			QMessageBox.information(None, 'Focal statistics', 'Choose input raster file!')
			return
		t1 = time()
		self.ui.lInfo.setText('Reading input file...')
		self.repaint()
		rLayer = getMapLayerByName(self.ui.cbLayers.currentText())
		if rLayer:
			path = str(rLayer.source())
		else:
			path = str(self.ui.cbLayers.currentText())
		
		gd = gdal.Open(path)
		if not gd:
			QMessageBox.information(None, 'Focal statistics', 'File not found!')
			return
			
		rasterInfo = focal.rasterInfo(gd)
		if self.ui.cbDataType.isChecked():
			rasterInfo.gdalDataType = focal.dataTypesFromInt[self.ui.cbType.currentIndex()]
		
		raster = gd.ReadAsArray()
		gd = None
		
		#QInputDialog.getText(None, '','', QLineEdit.Normal, str(rasterInfo.noDataValue))

		mask = self.createSimpleMask()
		
		self.ui.lInfo.setText('Calculating...')
		self.repaint()
		try:
			if self.ui.cbMethod.currentIndex() == 0:
				calcRast = focal.rasterSum(raster, mask, rasterInfo)
				raster = None
			elif self.ui.cbMethod.currentIndex() == 1:
				calcRast = focal.rasterAverage(raster, mask, rasterInfo)
			elif self.ui.cbMethod.currentIndex() == 2:
				calcRast = focal.rasterMedian(raster, mask, rasterInfo)
				raster = None
			elif self.ui.cbMethod.currentIndex() == 3:
				calcRast = focal.rasterMin(raster, mask, rasterInfo)
				raster = None
			elif self.ui.cbMethod.currentIndex() == 4:
				calcRast = focal.rasterMax(raster, mask, rasterInfo)
				raster = None
			elif self.ui.cbMethod.currentIndex() == 5:
				calcRast = focal.rasterRange(raster, mask, rasterInfo)
				raster = None
				
		except MemoryError:
			QMessageBox.critical(None, 'Focal Statistics', 'Memory Error')
			raster = None
			self.ui.lInfo.setText('Memory error...')
			self.repaint()
			return
		
		self.ui.lInfo.setText('Creating new file...')
		self.repaint()
		
		if self.ui.cbDataType.isChecked():
			focal.SaveToFile(str(self.ui.eOutput.text()), calcRast, rasterInfo, self.ui.cbType.currentIndex())
		else:
			focal.SaveToFile(str(self.ui.eOutput.text()), calcRast, rasterInfo)
			
		if self.ui.cbAddResult.isChecked():
			self.ui.lInfo.setText('Loading file...')
			self.repaint()
			r = QgsRasterLayer(str(self.ui.eOutput.text()), 'nowy_raster')
			QgsMapLayerRegistry.instance().addMapLayer(r)
		
		t2 = time()
		t = t2-t1
		self.ui.lInfo.setText('Finished')
		QMessageBox.information(None, 'Focal Statistics', 'Done in %.2f sec.' % t)


"""
def rasterRange(mask, rInfo, iterations=1):
	#noDataValue indexes
	global raster
	noDataArray = np.where(raster == rInfo.noDataValue)
	#calc MaxValue array
	rMax = focal.rasterMax(raster, mask, rInfo)
	
	#set noDataValue
	raster[noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	#calc MinValue array
	rMin = focal.rasterMin(raster, mask, rInfo)
	#print type(raster)#, raster
	raster = None #free some memory
	#del raster
	#gc.collect()
	#print type(raster)
	rRange = rMax-rMin #calc range
	rRange[noDataArray[0], noDataArray[1]] = rInfo.noDataValue #set noDataValue
	return rRange"""
	
	
		
def getLayersNames():
	layermap = QgsMapLayerRegistry.instance().mapLayers()
	layerlist = []
	#iteracja po elementach slownika
	for name, layer in layermap.iteritems(): 
		#sprawdzenie, czy dana warstwa jest warstwa rastrowa
		if layer.type() == QgsMapLayer.RasterLayer:
			#jesli tak, to dodaj do listy warstw
			layerlist.append( layer.name() )
	return layerlist #zwrocenie listy warstw
	
def getMapLayerByName(myName):
	#pobranie wszystkich wczytanych warstw
	layermap = QgsMapLayerRegistry.instance().mapLayers()
	for name, layer in layermap.iteritems():
		#sprawdzenie czy dana warstwa ma szukana nazwe
		if layer.name() == myName:
			return layer