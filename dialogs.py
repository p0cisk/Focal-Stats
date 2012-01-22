from PyQt4.QtGui import *
from PyQt4.QtCore import *

from os.path import splitext, dirname

def saveDialog(parent):
	"""Shows a save file dialog and return the selected file path."""
	settings = QSettings()
	key = '/UI/lastShapefileDir'
	outDir = settings.value(key).toString()

	filter = 'GeoTiff (*.tiff)'
	outFilePath = QFileDialog.getSaveFileName(parent, parent.tr('Save as GeoTiff'), outDir, filter)
	outFilePath = unicode(outFilePath)

	if outFilePath:
		root, ext = splitext(outFilePath)
		if ext.upper() != '.TIFF':
			outFilePath = '%s.tiff' % outFilePath
		outDir = dirname(outFilePath)
		settings.setValue(key, outDir)

	return outFilePath
	
def openDialog(parent):
	settings = QSettings()
	key = '/UI/lastShapefileDir'
	outDir = settings.value(key).toString()

	filter = 'GeoTiff (*.tif *.tiff)'
	outFilePath = QFileDialog.getOpenFileName(parent, parent.tr('Open raster file'), outDir, filter)
	outFilePath = unicode(outFilePath)
	
	if outFilePath:
		root, ext = splitext(outFilePath)
		outDir = dirname(outFilePath)
		settings.setValue(key, outDir)
	
	return outFilePath

def openDir(parent):
	settings = QSettings()
	key = '/UI/lastShapefileDir'
	outDir = settings.value(key).toString()

	outPath = QFileDialog.getExistingDirectory(parent, 'Focal Statistics', outDir)#, QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
	return outPath