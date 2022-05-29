# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from PyQt5.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterString,
                       QgsProcessingContext,
                       QgsProcessingFeedback,
                       QgsProcessingParameterBoolean,
                       QgsGraduatedSymbolRenderer,
                       QgsGradientColorRamp,
                       QgsProcessingUtils,
                       QgsVectorLayer
                       )
#from qgis.core import *
import qgis.utils
from osgeo import gdal, osr, ogr
import numpy as np
import glob
import gdal_merge as merge
import ogr2ogr
import os
import rsigdal
from PyQt5 import QtGui

class ExampleProcessingAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    wrs_path = 'wrs_path' ####
    wrs_row = 'wrs_row'    ###
    obs_pre = 'obs_pre'
    obs_post = 'obs_post'
    thres = 'thres'
    mask = "mask"
    count = "count"
    count_path = None
    
    def __init__(self):
        super().__init__()

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ExampleProcessingAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'myscript'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('２時期のLandsat-8データによる土地変状箇所抽出')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Example scripts')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'examplescripts'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        text = "このツールでは、２時期のLandsat-8の植生指数の差分から土地変状箇所を抽出する。"
        return self.tr(text)

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        
        self.addParameter(
            QgsProcessingParameterString(
                self.wrs_path,
                self.tr('WRS_Path (ex. 112)'),
                defaultValue = "112",
                optional = True
            )
        )
        
        self.addParameter(
            QgsProcessingParameterString(
                self.wrs_row,
                self.tr('WRS_Row (ex. 038)'),
                defaultValue = "038",
                optional = True
            )
        )
        
        self.addParameter(
            QgsProcessingParameterString(
                self.obs_pre,
                self.tr('1st observation: archive (ex. 20140502)'),
                defaultValue = 20140502,
                optional = True
            )
        )
        
        self.addParameter(
            QgsProcessingParameterString(
                self.obs_post,
                self.tr('2nd observation (ex. 20150521)'),
                defaultValue = 20150521,
                optional = True
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.thres,
                self.tr('thres'),
                type = QgsProcessingParameterNumber.Double,
                defaultValue = 0.3,
                optional = True
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.mask,
                self.tr('mask'),
                optional = True
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.count,
                self.tr('count'),
                optional = True
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        # Dummy function for thread safe mode
        def dummy(alg, context, feedback):
            pass

        tmp_dir = r"D:\tmp"
        sc_dir = r"D:\\Shortcuts"

        wrs_path = self.parameterAsString(parameters, self.wrs_path, context)
        wrs_row = self.parameterAsString(parameters, self.wrs_row, context)
        obs_pre = self.parameterAsString(parameters, self.obs_pre, context)
        obs_post = self.parameterAsString(parameters, self.obs_post, context)
        thres = self.parameterAsDouble(parameters, self.thres, context)
        mask_bool = self.parameterAsBool(parameters, self.mask, context)
        count_bool = self.parameterAsBool(parameters, self.count, context)
        self.count_bool = count_bool
        bandnumber = 1
        dtype = gdal.GDT_Float32

        proc = rsigdal.ProcessClass()

        path_list = proc.get_landsat8_result_path_list(wrs_path, wrs_row, obs_pre, obs_post)

        for dir_path in path_list:
            if not os.path.exists(dir_path):
                os.mkdir(dir_path)

        nir_pre_arr, red_pre_arr, mask_pre_arr = proc.read_landsat8_as_array(wrs_path, wrs_row, obs_pre)
        nir_post_arr, red_post_arr, mask_post_arr = proc.read_landsat8_as_array(wrs_path, wrs_row, obs_post)

        nir_pre_path = proc.get_landsat8_data_path(wrs_path, wrs_row, obs_pre, "band5")
        prop = proc.get_property_from_geotiff(nir_pre_path)

        mtl_path = proc.get_landsat8_txt_path(wrs_path, wrs_row, obs_pre, "MTL")
        meta_dict = proc.read_landsat8_mtl(mtl_path)
        ulx = meta_dict["CORNER_UL_PROJECTION_X_PRODUCT"]
        uly = meta_dict["CORNER_UL_PROJECTION_Y_PRODUCT"]
        lrx = meta_dict["CORNER_LR_PROJECTION_X_PRODUCT"]
        lry = meta_dict["CORNER_LR_PROJECTION_Y_PRODUCT"]
        baseEPSG = "EPSG:" + str(prop[4])

        maskfunc = proc.landsat8_mask_function(mask_pre_arr)
        nir_pre_arr[maskfunc] = np.nan
        red_pre_arr[maskfunc] = np.nan
        del mask_pre_arr

        ndvi_pre = proc.normalization(red_pre_arr, nir_pre_arr)
        proc.create_geotiff_from_array(ndvi_pre, os.path.join(tmp_dir, "PRE.tif"), bandnumber, dtype, np.nan, prop)

        del ndvi_pre, nir_pre_arr, red_pre_arr

        maskfunc = proc.landsat8_mask_function(mask_post_arr)
        nir_post_arr[maskfunc] = np.nan
        red_post_arr[maskfunc] = np.nan
        del mask_post_arr

        ndvi_post = proc.normalization(red_post_arr, nir_post_arr)

        nir_post_path = proc.get_landsat8_data_path(wrs_path, wrs_row, obs_post, "band5")
        prop = proc.get_property_from_geotiff(nir_post_path)
        proc.create_geotiff_from_array(ndvi_post, os.path.join(tmp_dir, "POST.tif"), bandnumber, dtype, np.nan, prop)

        del ndvi_post, nir_post_arr, red_post_arr

        if mask_bool:
            llon = float(min([meta_dict["CORNER_UL_LON_PRODUCT"], meta_dict["CORNER_LL_LON_PRODUCT"]]))
            hlon = float(max([meta_dict["CORNER_UR_LON_PRODUCT"], meta_dict["CORNER_LR_LON_PRODUCT"]]))
            llat = float(min([meta_dict["CORNER_LL_LAT_PRODUCT"], meta_dict["CORNER_LR_LAT_PRODUCT"]]))
            hlat = float(max([meta_dict["CORNER_UL_LAT_PRODUCT"], meta_dict["CORNER_UR_LAT_PRODUCT"]]))

            lc_path_list = proc.get_landcover_path(llon, hlon, llat, hlat)

            options = ['', '-o', os.path.join(tmp_dir, "landcover.tif")]

            for lc_path in lc_path_list:
                options.append(lc_path)
            merge.main(options)

            gdal.Warp(os.path.join(tmp_dir, "landcover2.tif"), os.path.join(tmp_dir, "landcover.tif"),
                      dstSRS=baseEPSG, xRes=15, yRes=15, resampleAlg="bilinear")

            options = ['', '-o', os.path.join(tmp_dir, "NDVI_merged.tif"),
                       '-separate', '-n', '0', '-ul_lr', str(ulx), str(uly), str(lrx), str(lry),
                       os.path.join(tmp_dir, "PRE.tif"), os.path.join(tmp_dir, "POST.tif"),
                       os.path.join(tmp_dir, "landcover2.tif")]

        else:
            options = ['', '-o', os.path.join(tmp_dir, "NDVI_merged.tif"),
                       '-separate', '-n', '0', '-ul_lr', str(ulx), str(uly), str(lrx), str(lry),
                       os.path.join(tmp_dir, "PRE.tif"), os.path.join(tmp_dir, "POST.tif")]

        if os.path.exists(os.path.join(tmp_dir, "NDVI_merged.tif")):
            os.remove(os.path.join(tmp_dir, "NDVI_merged.tif"))
        merge.main(options)

        ndvi_merged = gdal.Open(os.path.join(tmp_dir, "NDVI_merged.tif"), gdal.GA_ReadOnly)
        ndvi_pre_arr = ndvi_merged.GetRasterBand(1).ReadAsArray()
        ndvi_post_arr = ndvi_merged.GetRasterBand(2).ReadAsArray()

        ndvi_diff = ndvi_pre_arr - ndvi_post_arr
        ndvi_diff[np.isnan(ndvi_diff)] = 255

        if mask_bool:
            ndvi_mask_arr = ndvi_merged.GetRasterBand(3).ReadAsArray()
            ndvi_diff[(ndvi_diff != 255) & (ndvi_diff >= thres) & (ndvi_mask_arr >= 5)] = 1
            ndvi_diff[(ndvi_diff != 255) & (ndvi_diff < thres) & (ndvi_mask_arr >= 5)] = 0
        else:
            ndvi_diff[(ndvi_diff != 255) & (ndvi_diff >= thres)] = 1
            ndvi_diff[(ndvi_diff != 255) & (ndvi_diff < thres)] = 0

        merge.main(['', '-o', os.path.join(tmp_dir, "NDVI_merged.tif"),
                    '-separate', '-n', '0', os.path.join(tmp_dir, "PRE.tif"),
                    os.path.join(tmp_dir, "POST.tif")])


        prop = proc.get_property_from_geotiff(os.path.join(tmp_dir, "NDVI_merged.tif"))
        DIFF = proc.create_geotiff_from_array(ndvi_diff, os.path.join(path_list[3], "NDVI_diff.tif"), 1, gdal.GDT_Byte, 255, prop)

        ndvi_diff[ndvi_diff == 255] = 0
        DIFF_tmp = proc.create_geotiff_from_array(ndvi_diff, os.path.join(tmp_dir, "NDVI_diff.tif"), 1, gdal.GDT_Byte, 0, prop)

        del ndvi_merged, ndvi_pre_arr, ndvi_post_arr, ndvi_diff

        result_path = path_list[3] + os.sep + obs_pre + "_" + obs_post + ".shp"
        proc.raster_to_polygon(DIFF_tmp, os.path.splitext(result_path)[0])

        shapefile = ogr.Open(result_path, 1)
        layer = shapefile.GetLayer(0)
        new_field = ogr.FieldDefn("Area", ogr.OFTReal)
        new_field.SetWidth(20)
        new_field.SetPrecision(3)
        layer.CreateField(new_field)

        for feature in layer:
            geom = feature.GetGeometryRef()
            area = geom.GetArea()
            feature.SetField("Area", area)
            layer.SetFeature(feature)

        del shapefile, DIFF

        sc_file = obs_pre + "_" + obs_post + "_" + wrs_path + "_" + wrs_row + ".lnk"
        proc.create_shortcut(os.path.dirname(result_path), sc_dir, sc_file)

        if count_bool:
            polygon_dir = r"C:\RESTEC_QGIS_TOOLS\Share\region"
            polygon_path = glob.glob(polygon_dir + os.sep + "*.shp")
            result_ras_path = os.path.join(path_list[3], "NDVI_diff.tif")
            count_shp = obs_pre + "_" + obs_post + "_count.shp" 
            count_shp_path = os.path.join(path_list[3], count_shp)
            reproj_ras_path = os.path.join(path_list[3], "reproj_ras.tif")
            reproj_result_path = os.path.join(path_list[3], "reproj_result.shp")
	
            driver = gdal.GetDriverByName('GTIff')
            result_ras = gdal.Open(result_ras_path, gdal.GA_ReadOnly)
	
            shp = ogr.Open(polygon_path[0], 1)
            shp_lyr = shp.GetLayer(0)
            shp_srs = shp_lyr.GetSpatialRef()
            shp_wkt = shp_srs.ExportToWkt()
	
            gdal.Warp(reproj_ras_path, result_ras, dstSRS = shp_wkt)	
            reproj_ras = gdal.Open(reproj_ras_path, gdal.GA_ReadOnly)
            ext = proc.get_raster_extent(reproj_ras)
            del reproj_ras
		
            ogr2ogr.main(['', '-clipsrc', str(ext[0]), str(ext[1]), str(ext[2]), str(ext[3]), count_shp_path, polygon_path[0]])
	
            count_ds = ogr.Open(count_shp_path, 0)
            count_lyr = count_ds.GetLayer(0)
            count_srs = count_lyr.GetSpatialRef()
            count_ds.Destroy()
	
            proc.reproj_lyr(result_path, reproj_result_path, count_srs)
	
            proc.count_shp_in_poly(reproj_result_path, count_shp_path)

            self.count_path = count_shp_path

            context.addLayerToLoadOnCompletion(self.count_path,
                                           QgsProcessingContext.LayerDetails(name="count",
                                                                             project=context.project()))


            reproj_list = proc.get_file_list(path_list[3], "reproj*")
            for file in reproj_list:
                os.remove(file)

        context.addLayerToLoadOnCompletion(result_path,
                                           QgsProcessingContext.LayerDetails(name=obs_post,
                                                                             project=context.project()))

        return {}

    def postProcessAlgorithm(self, context, feedback):
#        output = context.getMapLayer(self.result)
        if self.count_bool:
            layer = QgsProcessingUtils.mapLayerFromString(self.count_path, context)
            print(layer)
            field = 'Count'
            color1 = QtGui.QColor('#FFFFFF')
            color2 = QtGui.QColor('#FF0000')

            renderer = QgsGraduatedSymbolRenderer()
            renderer.setClassAttribute(field)
            layer.setRenderer(renderer)
            layer.renderer().updateClasses(layer, QgsGraduatedSymbolRenderer.EqualInterval, 10)
            layer.renderer().updateColorRamp(QgsGradientColorRamp(color1, color2))
    #        path = r'C:\RESTEC_QGIS_TOOLS\Share\style\count.qml'
    #        layer.loadNamedStyle(path)
            layer.triggerRepaint()
    
        return {}