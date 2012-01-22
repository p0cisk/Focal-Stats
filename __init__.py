"""
/***************************************************************************
 FocalStats
                                 A QGIS plugin
 Raster statistics over a specified neighborhood
                             -------------------
        begin                : 2011-11-22
        copyright            : (C) 2011 by Piotr Pociask
        email                : NA
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""
def name():
    return "Focal Statistics"
def description():
    return "Raster statistics over a specified neighborhood"
def category(): 
    return "Raster"
def version():
    return "Version 0.1"
def icon():
    return "icon.png"
def qgisMinimumVersion():
    return "1.0"
def classFactory(iface):
    # load FocalStats class from file FocalStats
    from focalstats import FocalStats
    return FocalStats(iface)
