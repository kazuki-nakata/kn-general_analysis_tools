import glob
from osgeo import gdal
import numpy as np
import os
from ..helpers.misc import import_config

class LANDSAT8L2A():
    def __init__(self,zipfile):
        self.output=None
        
        self.params=import_config(config_path=r'conf/modis_params.yaml')
        
        basename = os.path.basename(zipfile)
        dirname = os.path.dirname(zipfile)
        self.product='L2A'

        bands = glob.glob(dirname + os.sep + basename + "\*_band?.tif")[0]
        for band in bands: 
            band_name=band #should change
            self.ds[band_name] = gdal.Open(band, gdal.GA_ReadOnly)
        mask_path = glob.glob(dirname + os.sep + basename + "\*_pixel_qa.tif")[0]
        self.ds["qa"] = gdal.Open(mask_path, gdal.GA_ReadOnly)
        
    def landsat8_mask_function(self,mask_arr):
        mask_arr=self.ds["qa"]
        maskfunc = np.where(
            (mask_arr != 322)
        )
        return maskfunc

    def landsat8_mask_function(self,mask_arr):
        mask_array=self.ds["qa"]
        maskfunc = np.where(
            (mask_arr == 324) |
            (mask_arr == 388) |
            (mask_arr == 836) |
            (mask_arr == 900) |
            (mask_arr == 1348) |
            (mask_arr == 328) |
            (mask_arr == 392) |
            (mask_arr == 840) |
            (mask_arr == 904) |
            (mask_arr == 1350) |
            (mask_arr == 336) |
            (mask_arr == 368) |
            (mask_arr == 400) |
            (mask_arr == 432) |
            (mask_arr == 848) |
            (mask_arr == 880) |
            (mask_arr == 912) |
            (mask_arr == 944) |
            (mask_arr == 1352) |
            (mask_arr == 352) |
            (mask_arr == 368) |
            (mask_arr == 416) |
            (mask_arr == 432) |
            (mask_arr == 480) |
            (mask_arr == 864) |
            (mask_arr == 880) |
            (mask_arr == 928) |
            (mask_arr == 944) |
            (mask_arr == 992)
        )
        return maskfunc

    def read_landsat8_mtl(self,path):
        mtl = open(path, "r")
        attribs = mtl.readlines()

        metadata = {}
        for ilines in attribs:
            eq = ilines.rstrip('\n').strip().split("=")
            if len(eq) == 2:
                metadata[eq[0].strip()] = eq[1].strip()

        return metadata