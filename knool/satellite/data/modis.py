from osgeo import gdal, osr
import numpy as np
import os
from ...helpers.misc import import_config
from ..algorithm import thermal_sensor_process
from ...geodata import geo_info, geo_io
from pyhdf.SD import SD, SDC


class MXD021KM:
    def __init__(self, hdfpath=None, hdfpath_mxd03=None, read=None):
        if hdfpath is None:
            return

        self.output = None

        self.params = import_config(config_path=os.path.dirname(__file__) + os.sep + "conf/modis_params.yaml")

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "Level 1B"
        self.ref_name = "EV_1KM_RefSB"
        self.em_name = "EV_1KM_Emissive"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        self.subdsProd = {}
        for i, val in enumerate(self.ds.GetSubDatasets()[0:10]):
            name = val[0].split(":")[4]
            self.subdsID[name] = i
            self.subdsProd[name] = "MXD021KM"

        if not os.path.exists(hdfpath_mxd03):
            print("Ancillary data in MXD021KM is loaded insted of MXD03 product")
            for i, val in enumerate(self.ds.GetSubDatasets()[10:17]):
                name = val[0].split(":")[4]
                self.subdsID[name] = i + 10
                self.subdsProd[name] = "MXD021KM"

        else:
            self.hdfpath_mxd03 = hdfpath_mxd03
            self.ds_mxd03 = gdal.Open(hdfpath_mxd03, gdal.GA_ReadOnly)
            for i, val in enumerate(self.ds_mxd03.GetSubDatasets()):
                name = val[0].split(":")[4]
                self.subdsID[name] = i
                self.subdsProd[name] = "MXD03"

        if read == "em":
            self.read_em()
        elif read == "ref":
            self.read_ref()

    def get_latlon_array(self):
        if self.hdfpath_mxd03 is not None:
            hdf_geo = SD(self.hdfpath_mxd03, SDC.READ)
            lat = hdf_geo.select("Latitude")
            latitude = lat[:, :]
            lon = hdf_geo.select("Longitude")
            longitude = lon[:, :]
            return latitude, longitude
        else:
            print("Specify MXD03 product together with MXD021KM")

    def get_subds(self, subdsID, product_name):
        if product_name == "MXD021KM":
            ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        elif product_name == "MXD03":
            ds = gdal.Open(self.ds_mxd03.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_ref_bnames(self):
        ds = gdal.Open(self.ds.GetSubDatasets()[self.subdsID[self.ref_name]][0], gdal.GA_ReadOnly)
        return ds.GetMetadata()["band_names"].split(",")

    def get_em_bnames(self):
        ds = gdal.Open(self.ds.GetSubDatasets()[self.subdsID[self.em_name]][0], gdal.GA_ReadOnly)
        return ds.GetMetadata()["band_names"].split(",")

    def read_em(self):
        self.ds_em = self.get_subds(self.subdsID[self.em_name], self.subdsProd[self.em_name])
        print("Succeeded. ds_em was incorporated inside the object")

    def read_ref(self):
        self.ds_ref = self.get_subds(self.subdsID[self.ref_name], self.subdsProd[self.ref_name])
        print("Succeeded. ds_ref was incorporated inside the object")

    def _radiance_to_bt(self, radiance, band_name):
        h = 6.626e-34  # planch constant
        c = 2.998e8  # velocity of light
        k = 1.3806e-23  # bolztmann constant
        lamda = float(self.params["wavelength"][int(band_name)]) * (1.0e-6)
        a = h * c / (k * lamda)
        b = 2 * h * c**2 / lamda**5
        print("coeffs -> band ", band_name, ":", a, b)
        return a / np.log(1 + b / radiance)

    def get_sensor_zenith(self):
        ds_sz = self.get_subds(self.subdsID["SensorZenith"], self.subdsProd["SensorZenith"])
        sensor_zenith = float(ds_sz.GetMetadata()["scale_factor"]) * ds_sz.ReadAsArray()
        return sensor_zenith

    def get_calibrated_em(self, band_min=0, band_max=10000):
        inArray = self.ds_em.ReadAsArray()
        if band_max == 10000:
            band_max = inArray.shape[0]

        offsets = np.array(self.ds_em.GetMetadata()["radiance_offsets"].split(",")).astype(np.float32)
        scales = np.array(self.ds_em.GetMetadata()["radiance_scales"].split(",")).astype(np.float32)
        print("offset_vals:", offsets[band_min:band_max])
        print("scale_vals:", scales[band_min:band_max])
        bt = np.array(
            [
                self._radiance_to_bt(
                    ((inArray[band, :, :] - offsets[band]) * scales[band]) * 1.0e6, self.get_em_bnames()[band]
                )
                for band in range(band_min, band_max)
            ]
        )
        self.output = bt
        self.output_prop = geo_info.get_property_from_raster_with_gcps(self.ds_em)
        return bt

    def get_calibrated_ref(self, band_min=0, band_max=10000):  # not test
        inArray = self.ds_ref.ReadAsArray()

        if band_min == 0 and band_max == 10000:
            band_min = 0
            band_max = inArray.shape[0]

        offsets = np.array(self.ds_ref.GetMetadata()["radiance_offsets"].split(",")).astype(np.float32)
        scales = np.array(self.ds_ref.GetMetadata()["radiance_scales"].split(",")).astype(np.float32)
        print("offset_vals:", offsets[band_min:band_max])
        print("scale_vals:", scales[band_min:band_max])
        ref = np.array([(inArray[band, :, :] - offsets[band]) * scales[band] for band in range(band_min, band_max)])
        self.output = ref
        self.output_prop = geo_info.get_property_from_raster_with_gcps(self.ds_ref)
        return ref

    def calc_IST(self, config_path=os.path.dirname(__file__) + os.sep + "conf/modis_ist.yaml"):
        cfg = import_config(config_path=config_path)
        north_lat = float(self.ds_em.GetMetadata()["NORTHBOUNDINGCOORDINATE"])

        if north_lat > 0:
            cfg2 = cfg["arctic"]
        else:
            cfg2 = cfg["antarctic"]

        calib = self.get_calibrated_em(band_min=10, band_max=12)
        scan_zenith = self.get_sensor_zenith()
        ist = thermal_sensor_process.calc_IST_by_split_window(calib[0, :, :], calib[1, :, :], scan_zenith, cfg2)
        self.output = ist
        self.output_prop = geo_info.get_property_from_raster_with_gcps(self.ds_em)
        return ist

    def set_output_prop(self, gcp_x=20, gcp_y=10):
        lat, lon = self.get_latlon_array()
        length, width = lat.shape

        gcps = []
        for ai in np.linspace(0, length - 1, gcp_y):
            for aj in np.linspace(0, width - 1, gcp_x):
                i = int(ai)
                j = int(aj)
                gcps.append(gdal.GCP(float(lon[i][j]), float(lat[i][j]), 0.0, j + 0.5, i + 0.5))

        source_ref = osr.SpatialReference()
        source_ref.ImportFromEPSG(4326)

        self.output_prop = [width, length, source_ref, gcps]

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")


class MXD29:
    def __init__(self, hdfpath=None):
        if hdfpath == None:
            return

        self.output = None

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "SeaIceProduct Level 2"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()[0:2]):
            name = val[0].split(":")[4]
            self.subdsID[name] = i

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_ist(self):
        ds_ist = self.get_subds(self.subdsID["Ice_Surface_Temperature"])
        scale_factor = float(ds_ist.GetMetadata()["scale_factor"])
        add_offset = float(ds_ist.GetMetadata()["add_offset"])
        ist_array = scale_factor * (ds_ist.ReadAsArray() - add_offset)
        self.output = ist_array
        self.output_prop = geo_info.get_property_from_raster_with_gcps(ds_ist)
        return ist_array

    def get_qa(self):
        ds_qa = self.get_subds(self.subdsID["Ice_Surface_Temperature_Pixel_QA"])
        self.output = ds_qa.ReadAsArray()
        self.output_prop = geo_info.get_property_from_raster_with_gcps(ds_qa)
        return self.output

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")


class MXD35_L2:
    def __init__(self, hdfpath=None):
        if hdfpath == None:
            return

        self.output = None

        self.basename = os.path.basename(hdfpath)
        self.dirname = os.path.dirname(hdfpath)

        self.product = "Cloud mask Level 2"

        self.ds = gdal.Open(hdfpath, gdal.GA_ReadOnly)
        self.subdsID = {}
        for i, val in enumerate(self.ds.GetSubDatasets()):
            name = val[0].split(":")[4]
            self.subdsID[name] = i

    def get_subds(self, subdsID):
        ds = gdal.Open(self.ds.GetSubDatasets()[subdsID][0], gdal.GA_ReadOnly)
        return ds

    def get_cmask(self):
        ds_cloud = self.get_subds(self.subdsID["Cloud_Mask"])
        cloud_array = ds_cloud.ReadAsArray()
        self.output = cloud_array
        self.output_prop = geo_info.get_property_from_raster_with_gcps(ds_cloud)
        return cloud_array

    def export_output(self, filepath, no_data=None, file_type="GTiff", dtype=gdal.GDT_Float32):
        geo_io.make_raster_with_gcps_from_array(self.output, filepath, dtype, no_data, self.output_prop, file_type)
        print("Exported")


class Open:
    def __new__(cls, hdfpath, *args):
        product_name = os.path.basename(hdfpath).split(".")[0]
        if product_name == "MOD021KM" or product_name == "MYD021KM":
            load_class = super().__new__(type(MXD021KM(None, None, None)))
            load_class.__init__(hdfpath, *args)
        elif product_name == "MOD35_L2" or product_name == "MYD35_L2":
            load_class = super().__new__(type(MXD35_L2(None)))
            load_class.__init__(hdfpath)
        elif product_name == "MOD29" or product_name == "MYD29":
            load_class = super().__new__(type(MXD29(None)))
            load_class.__init__(hdfpath)
        else:
            print("the product is not supported")
        return load_class
