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
"""
# Import the PyQt and QGIS libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from focalstatsdialog import FocalStatsDialog

import sys
sys.path.append('D:\eclipse\plugins\org.python.pydev.debug_2.2.4.2011110216\pysrc')
from pydevd import settrace

class FocalStats:

    def __init__(self, iface):
        # Save reference to the QGIS interface
        self.iface = iface

    def initGui(self):
        # Create action that will start plugin configuration
        self.action = QAction(QIcon(":/plugins/focalstats/icon.png"), \
                              "Focal Statistics", self.iface.mainWindow())
        self.actionAbout = QAction( QIcon( ":/about.png" ), "About",  \
                                    self.iface.mainWindow() )
        # connect the action to the run method
        QObject.connect(self.action, SIGNAL("triggered()"), self.run)
        QObject.connect( self.actionAbout, SIGNAL( "triggered()" ), self.about )

        # Add toolbar button and menu item
        
        if hasattr(self.iface, "addPluginToRasterMenu"): 
            # new menu available so add both actions into PluginName submenu 
            # under Raster menu 
            self.iface.addPluginToRasterMenu( "Focal Statistics", self.action )
            self.iface.addPluginToRasterMenu( "Focal Statistics", self.actionAbout )  
            # and add Run button to the Raster panel 
            self.iface.addRasterToolBarIcon( self.action )
             
        else: 
            # oops... old QGIS without additional menus. Place plugin under 
            # Plugins menu as usual 
            self.iface.addPluginToMenu( "Focal Statistics", self.action )
            self.iface.addPluginToMenu( "Focal Statistics", self.actionAbout )
            # and add Run button to the Plugins panel 
            self.iface.addToolBarIcon( self.action )
        
        #self.iface.addToolBarIcon(self.action)
        #self.iface.addPluginToMenu("&Focal Statistics", self.action)

    def unload(self):
        if hasattr(self.iface, "addPluginToRasterMenu"): 
            # new menu used, remove submenus from main Raster menu 
            self.iface.removePluginRasterMenu( "Focal Statistics", self.action )
            self.iface.removePluginRasterMenu( "Focal Statistics", self.actionAbout )  
            # also remove button from Raster toolbar 
            self.iface.removeRasterToolBarIcon( self.action ) 
        else: 
            # Plugins menu used, remove submenu and toolbar button 
            self.iface.removePluginMenu( "Focal Statistics", self.action )
            self.iface.removePluginMenu( "Focal Statistics", self.actionAbout )  
            self.iface.removeToolBarIcon( self.action )
            # Remove the plugin menu item and icon
        
        #self.iface.removePluginMenu("&Focal Statistics",self.action)
        #self.iface.removeToolBarIcon(self.action)

    # run method that performs all the real work
    def run(self):

        # create and show the dialog
        dlg = FocalStatsDialog()
        # show the dialog
        dlg.show()
        result = dlg.exec_()
        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code
            pass

    def about(self):
        QMessageBox.about(None, 'About...', 'Focal Statistics \n0.1 \n\nCreated by \nPiotr Pociask')
        