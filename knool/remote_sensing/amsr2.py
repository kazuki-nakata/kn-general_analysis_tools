from osgeo import gdal, osr
import numpy as np
import os
from ..helpers.misc import import_config
from ..geodata_processor import geo_info, geo_io


class AMSR2_L1R:
    def __init__(self, hdfpath=None):
        if hdfpath is None:
            return

        self.output = None

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "AMSR2_L1R"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[2][2:]
            self.subdsID[name] = i

        self.clip_array = slice(None)

    def keys(self):
        print(self.subdsID)

    def set_lat_range(self, lamin, lamax):
        lat = self.get_subds(self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
        lat_1d_min = np.min(lat, axis=1)
        lat_1d_max = np.max(lat, axis=1)
        self.clip_array = (lat_1d_min >= lamin) & (lat_1d_max <= lamax)

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_brightness_temperature(self, subdsID):
        return self.get_subds(subdsID).ReadAsArray()[self.clip_array] * 0.01

    def get_latlon(self, resolution="low"):
        # Positions of all resample data in L1R product are ajusted so as to they is consistent with latlons of 89A.
        lat = self.get_subds(self.subdsID["Latitude_of_Observation_Point_for_89A"]).ReadAsArray()
        lon = self.get_subds(self.subdsID["Longitude_of_Observation_Point_for_89A"]).ReadAsArray()
        if resolution == "low":
            lat = lat[:, 0::2]
            lon = lon[:, 0::2]
        lat = lat[self.clip_array]
        lon = lon[self.clip_array]
        return lat, lon

    def get_earth_azimuth(self):
        eaz = self.get_subds(self.subdsID["Earth_Azimuth"]).ReadAsArray()[self.clip_array] * 0.01
        self.output = eaz
        return eaz

    def get_land_ocean_flag(self, resolution="36GHz"):
        num = {"6GHz": 0, "10GHz": 1, "18GHz": 2, "36GHz": 3}
        loflag = self.get_subds(self.subdsID["Land_Ocean_Flag_6_to_36"]).ReadAsArray()[num[resolution]][self.clip_array]
        self.output = loflag
        return loflag

    def set_output_prop(self):
        lat, lon = self.get_latlon()
        length, width = lat.shape

        gcps = []
        for ai in np.linspace(0, length - 1, 5):
            for aj in np.linspace(0, width - 1, 10):
                i = int(ai)
                j = int(aj)
                # print(lat[i][0],lon[i][0], 0, i+1, 1)
                print(i, length, j, width, float(lon[i][j]), float(lat[i][j]))
                gcps.append(gdal.GCP(float(lon[i][j]), float(lat[i][j]), 0.0, j + 0.5, i + 0.5))

        # gcps.append(
        # gdal.GCP(float(lon[i][int(width / 2)]), float(lat[i][int(width / 2)]), 0, int(width / 2) + 0.5, i + 0.5)
        # )
        # gcps.append(gdal.GCP(float(lon[i][width - 1]), float(lat[i][width - 1]), 0, width - 1 + 0.5, i + 0.5))

        source_ref = osr.SpatialReference()
        source_ref.ImportFromEPSG(4326)

        self.output_prop = [width, length, source_ref, gcps]

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")


class Open:
    def __new__(cls, hdfpath, *args):
        product_name = os.path.basename(hdfpath).split("_")[0]
        if product_name == "GW1AM2":
            load_class = super().__new__(type(AMSR2_L1R(None)))
            load_class.__init__(hdfpath)
        else:
            print("the product is not supported")
        return load_class
