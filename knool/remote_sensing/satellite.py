import gdal
import glob
import numpy as np
import os
from read_modis import read_modis
from read_palsar import read_palsar
from read_landsat8 import read_landsat8

class read_satellite():
    def __new__( cls, infile, satellite):
        if satellite == 'modis':
            load_class=super().__new__( type( read_modis(hdfpath) ) )
            load_class.__init__(hdfpath)
        else:
            print('the satellite data is not supported')
        return load_class
    
    def import_output(self,in_array):
        self.output=in_array
        
    def export(self,output_path): #not test
        pass
