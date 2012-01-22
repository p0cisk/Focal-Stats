import numpy as np
from osgeo import gdal#, osr
from osgeo.gdalconst import *
#from time import time
#import cProfile
#import gc

import sys
sys.path.append('D:\eclipse\plugins\org.python.pydev.debug_2.2.4.2011110216\pysrc')
from pydevd import settrace

"""
funkcje:
-suma - ok
-srednia - 
-mediana - ok
-minimum - ok
-maksimum - ok
-zasieg - ok
-Majority ?
-Minority ?
-odchylenie standardowe
-wariacja

"""
"""
GDALDataType:
GDT_Unknown = 0, GDT_Byte = 1, GDT_UInt16 = 2, GDT_Int16 = 3, 
GDT_UInt32 = 4, GDT_Int32 = 5, GDT_Float32 = 6, GDT_Float64 = 7, 
GDT_CInt16 = 8, GDT_CInt32 = 9, GDT_CFloat32 = 10, GDT_CFloat64 = 11, 
GDT_TypeCount = 12
"""
	
dataTypes = {gdal.GDT_Byte : int, #@UndefinedVariable
			gdal.GDT_UInt16 : int, #@UndefinedVariable
			gdal.GDT_Int16 : int, #@UndefinedVariable
			gdal.GDT_UInt32 : int, #@UndefinedVariable
			gdal.GDT_Int32 : int, #@UndefinedVariable
			gdal.GDT_Float32 : float, #@UndefinedVariable
			gdal.GDT_Float64 : float, #@UndefinedVariable
			gdal.GDT_CInt16 : np.complex,
			gdal.GDT_CInt32 : np.complex,
			gdal.GDT_CFloat32 : np.complex,
			gdal.GDT_CFloat64 : np.complex
			}
			
dataTypesFromInt = {0 : gdal.GDT_Byte,
					1 : gdal.GDT_UInt16,
					2 : gdal.GDT_Int16,
					3 : gdal.GDT_UInt32,
					4 : gdal.GDT_Int32,
					5 : gdal.GDT_Float32,
					6 : gdal.GDT_Float64,
					7 : gdal.GDT_CInt16,
					8 : gdal.GDT_CInt32,
					9 : gdal.GDT_CFloat32,
					10 : gdal.GDT_CFloat64
					}

class rasterInfo(object):
	def __init__(self, gd=None, bandNr=1):
		if gd:
			band = gd.GetRasterBand(bandNr)
			
			self.maxValue = band.GetMaximum()
			self.minValue = band.GetMinimum()

			
			if self.maxValue is None or self.minValue is None:
				(self.minValue, self.maxValue) = band.ComputeRasterMinMax(1)
				
			self.meanValue = (float(self.maxValue)+self.minValue)/2
			
			self.noDataValue = band.GetNoDataValue()
			if self.noDataValue is None: self.noDataValue = self.minValue-1
			
			self.dataType = dataTypes[band.DataType]
			self.gdalDataType = band.DataType
			
			self.projection = gd.GetProjection()
			self.geoTransform = gd.GetGeoTransform()
			
			
def rasterSum(raster, mask, rInfo, type=None, iterations=1):
	u = raster
	noDataArray = np.where(u==rInfo.noDataValue)
	u[u==rInfo.noDataValue]=0
	
		
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	wm = dx*2
	hm = dy*2
	
	
	for x in xrange(dx):
		u = np.insert(u, 0, 0, axis=1) #@UndefinedVariable
		u = np.insert(u, u.shape[1], 0, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, 0, axis=0)
		u = np.insert(u, u.shape[0], 0, axis=0)
	
	if type:
		out = np.zeros(u.shape, dtype=type)
	else:
		out = np.zeros(u.shape, dtype=rInfo.dataType)
		
	for i in xrange(iterations): #przy iteracji dodawane sa caly czas te same wartosci:/
		for x in xrange(wm):
			for y in xrange(hm):
				if mask[y,x] != 0:
					out[dy:-dy, dx:-dx] += (u[y:y-(dy+1*dy), x:x-(dx+1*dx)])#*np.array([mask[y,x]])

		for x in xrange(wm):
			if mask[-1,x] != 0:
				out[dy:-dy, dx:-dx] += u[(dy+1*dy):, x:x-(dx+1*dx)]#*np.array([mask[-1,x]])

		for y in xrange(hm):
			if mask[y,-1] != 0:
				out[dy:-dy, dx:-dx] += u[y:y-(dy+1*dy), (dx+1*dx):]#*np.array([mask[y,-1]])

		if mask[-1,-1] != 0:
			out[dy:-dy, dx:-dx] += u[(dy+1*dy):, (dx+1*dx):]#*np.array([mask[-1,-1]])

	print out
	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	
	return out[dy:-dy,dx:-dx]
	
	
def rasterMax(raster, mask, rInfo, iterations=1):
	u = raster
	noDataArray = np.where(u == rInfo.noDataValue)
	u[u==rInfo.noDataValue] = rInfo.minValue
	
	
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	#wm = dx*2
	#hm = dy*2
	
	for x in xrange(dx):
		u = np.insert(u, 0, rInfo.minValue, axis=1)
		u = np.insert(u, u.shape[1], rInfo.minValue, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, rInfo.minValue, axis=0)
		u = np.insert(u, u.shape[0], rInfo.minValue, axis=0)


	out = np.zeros(u.shape, dtype=rInfo.dataType)
	out[:,:] = rInfo.minValue
	
	for x in xrange(mask.shape[1]):
		for y in xrange(mask.shape[0]):
			if mask[-(y+1),-(x+1)] != 0:
				out[dy:-dy, dx:-dx] = np.maximum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])

	
	

	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	rOut = out[dy:-dy,dx:-dx]
	del out
	return rOut
	
	
def rasterMin(raster, mask, rInfo, iterations=1):
	u = raster
	noDataArray = np.where(u == rInfo.noDataValue)
	u[u==rInfo.noDataValue] = rInfo.maxValue
	
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	#wm = dx*2
	#hm = dy*2
	#print u
	
	for x in xrange(dx):
		u = np.insert(u, 0, rInfo.maxValue, axis=1)
		u = np.insert(u, u.shape[1], rInfo.maxValue, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, rInfo.maxValue, axis=0)
		u = np.insert(u, u.shape[0], rInfo.maxValue, axis=0)

	out = np.zeros(u.shape, dtype=rInfo.dataType)
	out[:,:] = rInfo.maxValue


	for x in xrange(mask.shape[1]):
		for y in xrange(mask.shape[0]):
			if mask[-(y+1),-(x+1)] != 0:
				out[dy:-dy, dx:-dx] = np.minimum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])

	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	rOut = out[dy:-dy,dx:-dx]
	del out
	return rOut

	
	
def rasterRange2(raster, mask, rInfo, iterations=1):
	u = raster
	noDataArray = np.where(u == rInfo.noDataValue)
	#u[u==rInfo.noDataValue] = rInfo.meanValue
	
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	#wm = dx*2
	#hm = dy*2
	#print u
	
	for x in xrange(dx):
		u = np.insert(u, 0, rInfo.noDataValue, axis=1)
		u = np.insert(u, u.shape[1], rInfo.noDataValue, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, rInfo.noDataValue, axis=0)
		u = np.insert(u, u.shape[0], rInfo.noDataValue, axis=0)

	out = np.zeros(u.shape, dtype=rInfo.dataType)
	out[:,:] = rInfo.noDataValue
	
	for x in xrange(mask.shape[1]):
		for y in xrange(mask.shape[0]):
			if mask[-(y+1),-(x+1)] != 0:
				out[dy:-dy, dx:-dx] = np.maximum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])-np.minimum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])


	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	return out[dy:-dy,dx:-dx]
	
	
def rasterRange(raster, mask, rInfo, iterations=1):
	#noDataValue indexes
	noDataArray = np.where(raster == rInfo.noDataValue)
	#calc MaxValue array
	rMax = rasterMax(raster, mask, rInfo)
	#print rInfo
	
	#set noDataValue
	raster[noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	#calc MinValue array
	rMin = rasterMin(raster, mask, rInfo)
	#print type(raster)#, raster
	
	rRange = rMax-rMin #calc range

	rRange[noDataArray[0], noDataArray[1]] = rInfo.noDataValue #set noDataValue
	return rRange
	
	
def rasterMedian(raster, mask, rInfo, iterations=1):
	u = raster
	noDataArray = np.where(u==rInfo.noDataValue)
	u[u==rInfo.noDataValue]=rInfo.noDataValue
		
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	#wm = dx*2
	#hm = dy*2
	
	for x in xrange(dx):
		u = np.insert(u, 0, rInfo.noDataValue, axis=1)
		u = np.insert(u, u.shape[1], rInfo.noDataValue, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, rInfo.noDataValue, axis=0)
		u = np.insert(u, u.shape[0], rInfo.noDataValue, axis=0)
		
	(uh, uw) = u.shape
	#print u, hm, wm
	
	out = np.zeros(u.shape, dtype=rInfo.dataType)
		
	for w in xrange(uw-dx-1):
		for h in xrange(uh-dy-1):
			s = u[h:h+mh,w:w+mw]#*mask
			s = s[np.where(s!=rInfo.noDataValue)]
			if len(s) == 0: continue #if all values are noDataValue then continue with next cell
			
			p = list(s)
			p.sort()
			l = len(p)
			l2 = l/2

			if l%2:
				m = p[l2]
			else:
				m = (float(p[l2])+p[(l2)-1])/2
			out[h+dy,w+dx] = m
		
		
	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	
	return out[dy:-dy,dx:-dx]
	
	
def rasterMedian2(raster, mask, rInfo, iterations=1):
	u = raster
	noDataArray = np.where(u==rInfo.noDataValue)
	u[u==rInfo.noDataValue]=rInfo.noDataValue
		
	(mh, mw) = mask.shape #wymiary maski uzytkownika

	dx = mw//2
	dy = mh//2
	#wm = dx*2
	#hm = dy*2
	
	for x in xrange(dx):
		u = np.insert(u, 0, rInfo.noDataValue, axis=1)
		u = np.insert(u, u.shape[1], rInfo.noDataValue, axis=1)

	for y in xrange(dy):
		u = np.insert(u, 0, rInfo.noDataValue, axis=0)
		u = np.insert(u, u.shape[0], rInfo.noDataValue, axis=0)
		
	(uh, uw) = u.shape
	#print u, hm, wm
	
	out = np.zeros(u.shape, dtype=rInfo.dataType)
		
	"""for w in xrange(uw-dx-1):
		for h in xrange(uh-dy-1):
			s = u[h:h+mh,w:w+mw]#*mask
			s = s[np.where(s!=rInfo.noDataValue)]
			if len(s) == 0: continue #if all values are noDataValue then continue with next cell
			
			p = list(s)
			p.sort()
			l = len(p)
			l2 = l/2

			if l%2:
				m = p[l2]
			else:
				m = (float(p[l2])+p[(l2)-1])/2
			out[h+dy,w+dx] = m"""
	out = focal_median.median_calc(u, uw, dx, uh, dy, mh, mw, out, rInfo.noDataValue)

		
		
	out[dy:-dy,dx:-dx][noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	
	return out[dy:-dy,dx:-dx]
	
	
def rasterAverage(raster, mask, rInfo, iterations=1):
	noDataArray = np.where(raster == rInfo.noDataValue)
	count = simpleMaskCount(raster.shape, mask.shape, noDataArray)

	out = rasterSum(raster, mask, rInfo)/count
	
	out[noDataArray[0], noDataArray[1]] = rInfo.noDataValue
	return out 


def simpleMaskCount(rasterShape, maskShape, noDataArray):
	(rY, rX) = rasterShape
	#rY = rasterShape[0]
	(mY, mX) = maskShape
	#mY = maskShape[0]
	
	out = np.zeros(rY*rX, dtype=int).reshape(rY,rX)
	count = mY*mX
	
	dx = mX/2
	dy = mY/2
	#wm = mX-1
	#hm = mY-1
	#print dx, dy, wm, hm
	
	out[dy:-dy,dx:-dx] = count
	

	for x in xrange(dx):
		out[dy:-dy,x] = mY*(dx+1+x)
		out[dy:-dy,-x-1] = mY*(dx+1+x)
		

	for y in xrange(dy):
		out[y,dx:-dx] = (dy+1+y)*mX
		out[-y-1,dx:-dx] = (dy+1+y)*mX
		
	
	for x in xrange(dx):
		for y in xrange(dy):
			val = (dy+1+y)*(dx+1+x)
			out[y, x] = val
			out[-(y+1),-(x+1)] = val
			out[y, -(x+1)] = val
			out[-(y+1),x] = val

	tmp = np.zeros((rY+dy*2)*(rX+dx*2), dtype=float).reshape(rY+dy*2,rX+dx*2)
	tmp[dy:-dy,dx:-dx] = out
	
	for i in xrange(len(noDataArray[0])):
		tmp[noDataArray[0][i]:noDataArray[0][i]+(2*dy)+1,noDataArray[1][i]:noDataArray[1][i]+(2*dx)+1] -= 1
		
	
	
	"""for x in range(dx):
		c = list(self.mask[:, dx-x:].flatten()).count(0)
		out[hm:-hm,x+dx] = self.mask[:, dx-x:].size-c
		#print c
		
		c = list(self.mask[:, :dx+1+x].flatten()).count(0)
		out[hm:-hm,-(x+1+dx)] = self.mask[:, :dx+1+x].size-c"""
	
	return tmp[dy:-dy,dx:-dx]
	

def mask_count(self):
	u = self.big_r#np.array(self.big_r, copy=True)
		
	dx = self.size_x
	dy = self.size_y
	wm = dx*2
	hm = dy*2
	#print u
	out = np.zeros(u.shape, dtype=float) #zmienic na float

	c = list(self.mask.flatten()).count(0)
	out[hm:-hm,wm:-wm] = self.mask.size-c


	#boki
	for x in range(dx):
		c = list(self.mask[:, dx-x:].flatten()).count(0)
		out[hm:-hm,x+dx] = self.mask[:, dx-x:].size-c
		#print c
		
		c = list(self.mask[:, :dx+1+x].flatten()).count(0)
		out[hm:-hm,-(x+1+dx)] = self.mask[:, :dx+1+x].size-c
		#gora, dol
	for y in range(dy):
		c = list(self.mask[dy-y:, :].flatten()).count(0)
		out[y+dy,wm:-wm] = self.mask[dy-y:, :].size-c

		c = list(self.mask[:dy+1+y, :].flatten()).count(0)
		out[-(y+1+dy),wm:-wm] = self.mask[:dy+1+y, :].size-c
		#rogi
	for x in range(dx):
		for y in range(dy):
			#lewy gorny
			c = list(self.mask[dy:, dx-x:].flatten()).count(0)
			out[y+dy, x+dx] = self.mask[dy:, dx-x:].size-c

			#prawy dolny
			c = list(self.mask[:dy+y+1, :dx+x+1].flatten()).count(0)
			out[-(y+1+dy),-(x+1+dx)] = self.mask[:dy+y+1, :dx+x+1].size-c
			#lewy dolny
			c = list(self.mask[:dy+y+1, dx-x:].flatten()).count(0)
			out[-(y+1+dy), x+dx] = self.mask[:dy+y+1, dx-x:].size-c

			#prawy gorny
			c = list(self.mask[dy:, :dx+x+1].flatten()).count(0)
			out[y+dy, -(x+1+dx)] = self.mask[dy:, :dx+x+1].size-c

		
	rast_nan = np.where(u[dy:-dy,dx:-dx]==self.noDataValue)
	rast_nan = (rast_nan[0]+dy,rast_nan[1]+dx) #dopasowanie do duzego rastra
	rev_mask = np.fliplr(self.mask[::-1])


	for i in range(rast_nan[0].size):
		ry = rast_nan[0][i]
		rx = rast_nan[1][i]
		out[ry-dy:ry+dy+1,rx-dx:rx+dx+1] -= rev_mask


	return out[dy:-dy,dx:-dx]
	
			
class focal(object):
	def __init__(self, mask, raster, rInfo):
		self.noDataValue = rInfo.noDataValue
		self.dataType = rInfo.dataType
		self.maxValue = rInfo.maxValue
		self.minValue = rInfo.minValue

		
		#print self.noDataValue, self.dataType, self.maxValue, self.minValue
		
		self.mask = mask #user mask
		(self.mh, self.mw) = self.mask.shape #wymiary maski uzytkownika

		self.size_x = self.mw//2
		self.size_y = self.mh//2

		tmp = np.array(raster, dtype=self.dataType)
		#dodanie dodatkowych osi zgodnie z rozmiarem maski uzytkownia i przygotowanie do stworzenia mask array
		for x in xrange(self.size_x):
			tmp = np.insert(tmp, 0, self.noDataValue, axis=1)
			tmp = np.insert(tmp, tmp.shape[1], self.noDataValue, axis=1)

		for y in xrange(self.size_y):
			tmp = np.insert(tmp, 0, self.noDataValue, axis=0)
			tmp = np.insert(tmp, tmp.shape[0], self.noDataValue, axis=0)
		#print tmp
		self.big_r = tmp #raster z dodatkowymi osiami
		#print self.big_r
		(self.h, self.w) = self.big_r.shape #nowy rozmiar


	def sum(self):
		u = self.big_r#np.array(self.big_r, copy=True)
		u[u==self.noDataValue]=0
		
		mask = self.mask
		dx = self.size_x
		dy = self.size_y
		wm = dx*2
		hm = dy*2

		out = np.zeros(u.shape, dtype=self.dataType)
		#print u

		for x in range(wm):
			for y in range(hm):
				if mask[y,x] != 0:
					out[dy:-dy, dx:-dx] += (u[y:y-(dy+1*dy), x:x-(dx+1*dx)])#*np.array([mask[y,x]])

		for x in range(wm):
			if mask[-1,x] != 0:
				out[dy:-dy, dx:-dx] += u[(dy+1*dy):, x:x-(dx+1*dx)]#*np.array([mask[-1,x]])

		for y in range(hm):
			if mask[y,-1] != 0:
				out[dy:-dy, dx:-dx] += u[y:y-(dy+1*dy), (dx+1*dx):]#*np.array([mask[y,-1]])

		if mask[-1,-1] != 0:
			out[dy:-dy, dx:-dx] += u[(dy+1*dy):, (dx+1*dx):]#*np.array([mask[-1,-1]])

		#print out[dy:-dy,dx:-dx]
		#print out
		out[self.big_r==self.noDataValue] = self.noDataValue
		#print out
		#print self.big_r
		return out[dy:-dy,dx:-dx]

	def average(self):
		out = np.array(self.sum()/self.mask_count(), dtype=float)
		out[self.big_r[self.size_y:-self.size_y, self.size_x:-self.size_x]==self.noDataValue] = self.noDataValue
		return out

	def max(self):
		u = self.big_r#np.array(self.big_r, copy=True)
		u[u==self.noDataValue]=self.minValue
		#print u
		
		mask = self.mask
		dx = self.size_x
		dy = self.size_y
		#wm = dx*2
		#hm = dy*2

		out = np.zeros(u.shape, dtype=self.dataType)
		for x in range(mask.shape[1]):
			for y in range(mask.shape[0]):
				if mask[-(y+1),-(x+1)] != 0:
					out[dy:-dy, dx:-dx] = np.maximum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])

		out[self.big_r==self.noDataValue] = self.noDataValue
		return out[dy:-dy,dx:-dx]

	def min(self):
		u = self.big_r#np.array(self.big_r, copy=True)
		u[u==self.noDataValue]=self.maxValue
		#print u
		
		mask = self.mask
		dx = self.size_x
		dy = self.size_y
		#wm = dx*2
		#hm = dy*2

		out = np.zeros(u.shape, dtype=self.dataType)
		out[:,:] = self.maxValue#np.nan_to_num(u).max()+1
		#print u
		#u[np.isnan(u)] = np.inf
		#print out
		for x in range(self.size_x):
			u[:,x] = u.max()+1
			u[:,-(1+x)] = u.max()+1
		for y in range(self.size_y):
			u[y,:] = u.max()+1
			u[-(1+y),:] = u.max()+1


		for x in range(mask.shape[1]):
			for y in range(mask.shape[0]):
				if mask[-(y+1),-(x+1)] != 0:
					out[dy:-dy, dx:-dx] = np.minimum(out[dy:-dy, dx:-dx], u[dy+(dy-y):u.shape[0]-y, dx+(dx-x):u.shape[1]-x])

		out[self.big_r==self.noDataValue] = self.noDataValue
		return out[dy:-dy,dx:-dx]

	def raster_range(self):
		return self.max()-self.min()

	def mask_count(self):
		u = self.big_r#np.array(self.big_r, copy=True)
		#u[u==self.noDataValue]=0
		
		dx = self.size_x
		dy = self.size_y
		wm = dx*2
		hm = dy*2
		
		#print u

		out = np.zeros(u.shape, dtype=float) #zmienic na float
		
		#print out
		

		#srodek
		#print self.mask
		c = list(self.mask.flatten()).count(0)
		out[hm:-hm,wm:-wm] = self.mask.size-c
		

		#boki
		for x in range(dx):
			c = list(self.mask[:, dx-x:].flatten()).count(0)
			out[hm:-hm,x+dx] = self.mask[:, dx-x:].size-c
			#print c
			
			c = list(self.mask[:, :dx+1+x].flatten()).count(0)
			out[hm:-hm,-(x+1+dx)] = self.mask[:, :dx+1+x].size-c

		#gora, dol
		for y in range(dy):
			c = list(self.mask[dy-y:, :].flatten()).count(0)
			out[y+dy,wm:-wm] = self.mask[dy-y:, :].size-c

			c = list(self.mask[:dy+1+y, :].flatten()).count(0)
			out[-(y+1+dy),wm:-wm] = self.mask[:dy+1+y, :].size-c

		#rogi
		for x in range(dx):
			for y in range(dy):
				#lewy gorny
				c = list(self.mask[dy:, dx-x:].flatten()).count(0)
				out[y+dy, x+dx] = self.mask[dy:, dx-x:].size-c

				#prawy dolny
				c = list(self.mask[:dy+y+1, :dx+x+1].flatten()).count(0)
				out[-(y+1+dy),-(x+1+dx)] = self.mask[:dy+y+1, :dx+x+1].size-c

				#lewy dolny
				c = list(self.mask[:dy+y+1, dx-x:].flatten()).count(0)
				out[-(y+1+dy), x+dx] = self.mask[:dy+y+1, dx-x:].size-c

				#prawy gorny
				c = list(self.mask[dy:, :dx+x+1].flatten()).count(0)
				out[y+dy, -(x+1+dx)] = self.mask[dy:, :dx+x+1].size-c

		#print out
		#if list(self.mask.flatten()).count(0) != 0:
		#rast_nan = np.where(np.isnan(u[dy:-dy,dx:-dx]))
		#rast_nan = (rast_nan[0]+dy,rast_nan[1]+dx) #dopasowanie do duzego rastra
		#rev_mask = np.fliplr(mask[::-1])
			
		rast_nan = np.where(u[dy:-dy,dx:-dx]==self.noDataValue)
		#print rast_nan
		rast_nan = (rast_nan[0]+dy,rast_nan[1]+dx) #dopasowanie do duzego rastra
		rev_mask = np.fliplr(self.mask[::-1])


			
		for i in range(rast_nan[0].size):
			ry = rast_nan[0][i]
			rx = rast_nan[1][i]
			out[ry-dy:ry+dy+1,rx-dx:rx+dx+1] -= rev_mask

		#out[self.big_r==self.noDataValue] = 1
		#print out[dy:-dy,dx:-dx]
		return out[dy:-dy,dx:-dx]
		
	
	def median1(self): #przerobic na maski (np.ma.median(mask)) ma.masked_invalid
		u = self.big_r
		#mask = self.mask
		
		dx = self.size_x
		dy = self.size_y
		wm = (dx*2)+1
		hm = (dy*2)+1
		
		out = np.zeros(u.shape, dtype=float)
		
		for w in range(self.w-dx-1):
			for h in range(self.h-dy-1):
				s = u[h:h+hm,w:w+wm]#*mask
				s = s[np.where(s>0)]
				out[h+dy,w+dx] = np.median(s)
		
		
		return np.array(out[dy:-dy,dx:-dx], dtype=float)
		
		
	def median2(self): #przerobic na maski (np.ma.median(mask)) ma.masked_invalid - wolne!!!
		u = self.big_r
		#mask = self.mask
		dx = self.size_x
		dy = self.size_y
		wm = (dx*2)+1
		hm = (dy*2)+1
		
		#out = np.zeros(u.shape, dtype=float)
		l = []
		#print u
		for x in range(mask.shape[1]-1):
			for y in range(mask.shape[0]-1):
				#print y, x
				l.append(u[y:-(hm-y)+1,x:-(wm-x)+1])
		
		for y in range(mask.shape[0]-1):
			l.append(u[y:-(hm-1-y),wm-1:])
			
		for x in range(mask.shape[1]-1):
			l.append(u[hm-1:,x:-(wm-1-x)])		
			
		l.append(u[hm-1:,wm-1:])
		
		#print l
		out = np.array(l, dtype=float)
		mout = np.ma.masked_invalid(out)
		
		"""
		for w in range(self.w-dx-1):
			for h in range(self.h-dy-1):
				s = u[h:h+hm,w:w+wm]#*mask
				s = s[np.where(s>0)]
				out[h+dy,w+dx] = np.median(s)"""
		
		
		return np.ma.median(mout)#np.array(out[dy:-dy,dx:-dx], dtype=float)
		
	def median3(self): #prawie 2x szybsze od median1
		u = self.big_r
		#mask = self.mask
		dx = self.size_x
		dy = self.size_y
		wm = (dx*2)+1
		hm = (dy*2)+1
		
		out = np.zeros(u.shape, dtype=float)
		
		for w in range(self.w-dx-1):
			for h in range(self.h-dy-1):
				s = u[h:h+hm,w:w+wm]#*mask
				s = s[np.where(s>0)]

				p = list(s)
				p.sort()
				l = len(p)
				l2 = l/2

				if l%2:
					m = p[l2]
				else:
					m = (float(p[l2])+p[(l2)-1])/2
				out[h+dy,w+dx] = m
		
		
		return np.array(out[dy:-dy,dx:-dx], dtype=float)



def SaveToFile(FileName, raster, rInfo, type=None): # do zrobienia
	format = "GTiff"
	driver = gdal.GetDriverByName( format )
	
	#print raster

	if type != None:
		dst_ds = driver.Create( FileName, raster.shape[1], raster.shape[0], 1, dataTypesFromInt[type], ['COMPRESS=DEFLATE'] )
	else:
		dst_ds = driver.Create( FileName, raster.shape[1], raster.shape[0], 1, rInfo.gdalDataType, ['COMPRESS=DEFLATE'] )
	
	
	dst_ds.SetGeoTransform( rInfo.geoTransform )
	
	dst_ds.SetProjection( rInfo.projection )
	
	#noDataValue = gd.GetRasterBand(1).GetNoDataValue()
	if rInfo.noDataValue: dst_ds.GetRasterBand(1).SetNoDataValue( rInfo.noDataValue )

	dst_ds.GetRasterBand(1).WriteArray( raster )

	# Once we're done, close properly the dataset
	dst_ds = None


if __name__ == '__main__':
	#liczba kolumn i wierszy musi byc nieparzysta!!!
	mask = np.array([[1, 1, 1],
					[1, 1, 1],
					[1, 1, 1]])

	#rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/gis-support/focal/dem.tif'
	#rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/gis-support/focal/model.tif'
	#rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/gis-support/lineamenty/MODEL_21.tif'
	#rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/swinka/dem.tif'
	#rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/cut_dem.asc'
	rasterFile = 'C:/Documents and Settings/Pocisk/Pulpit/stmp.tif'
	gd = gdal.Open(rasterFile)
	rInfo = rasterInfo(gd)
	#print dir(gd)
	#print gd.GetRasterBand(1).GetNoDataValue()
	
	raster = gd.ReadAsArray() #array with raster
	
	#r = rasterAverage(raster, mask, rInfo)
	#print r
	
	gd = None
	#nd = gd.GetRasterBand(1).GetNoDataValue()
	#dt = gd.GetRasterBand(1).DataType #dataType

	#print raster
	#raster[raster==gd.GetRasterBand(1).GetNoDataValue()]=np.nan
	#print raster
	#print np.where(np.isnan(raster))

	#f = focal(mask, raster, rInfo)
	#print 'sum'
	#cProfile.run('f.median()', 'C:/Documents and Settings/Pocisk/Pulpit/median.profile')
	#cProfile.run('f.median3()', 'C:/Documents and Settings/Pocisk/Pulpit/median3.profile')
	a = rasterSum(raster, mask, rInfo)
	#print raster
	#print 'Wynik:'
	#print a
	SaveToFile('C:/Documents and Settings/Pocisk/Pulpit/plik1.tif', a, rInfo)
	print 'ok'