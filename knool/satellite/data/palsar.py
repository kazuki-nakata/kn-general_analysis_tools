from osgeo import gdal
import numpy as np
import zipfile
import os
import pandas as pd
import xarray as xr

class PALSARL15():    
    def __init__(self,zip_file):
        self.product='l15'
        with zipfile.ZipFile(zip_file, 'r') as zip_data:
            infos = zip_data.infolist()
            for i,info in enumerate(infos):
                if info.filename.find("IMG") != -1:
                    img = info.filename
                elif(info.filename.find("/workreport")) != -1 :
                    meta_file = info.filename
            path = "/vsizip" + os.sep + zip_file + os.sep + img
            self.ds=gdal.Open(path,gdal.GA_ReadOnly)
            
            ds_name=["name","value"]
            with zip_data.open(meta_file) as f:
                self.df = pd.read_csv(f,delimiter="=",names=ds_name)

    def to_db(self):
        offset=-83.0
        return np.where(array!=0, 10*np.log10(np.square(array.astype(np.float32)))+offset, np.NaN)
    
class read_palsar():
    def __new__( cls,zip_file ):
        product_name=os.path.basename(hdfpath).split('.')[0]
        if product_name == 'MOD021KM' or product_name == 'MYD021KM':
            load_class=super().__new__( type( PALSARL15(hdfpath) ) )
            load_class.__init__(hdfpath)
        else:
            print('the product is not supported')
        return load_class