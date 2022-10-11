from osgeo import gdal, osr, ogr
import numpy as np
import os

OGRTypes = {int: ogr.OFTInteger, str: ogr.OFTString, float: ogr.OFTReal}


def copy_dataset(ds, outfile="/vsimem/output.tif"):
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.CreateCopy(outfile, ds)

    data = ds.ReadAsArray()

    if len(data.shape) == 2:
        num_band = 1
    else:
        num_band = data.shape[0]

    if num_band == 1:
        src_band = ds.GetRasterBand(1)
        band = dst_ds.GetRasterBand(1)
        band.WriteArray(data)
        band.SetNoDataValue(src_band.GetNoDataValue())
        band.FlushCache()
    else:
        for i in range(0, num_band):
            src_band = ds.GetRasterBand(i + 1)
            band = dst_ds.GetRasterBand(i + 1)
            band.WriteArray(data[i, :, :])
            band.SetNoDataValue(src_band.GetNoDataValue())
            band.FlushCache()
    return dst_ds


# fmt: off
def make_raster_from_array(data, filepath, pixel_x, pixel_y, num_band, dtype, no_data, file_type, geomode="proj",
                           geotrans=None, geoproj=None, gcps=None, gcpsrc=None):
    # fmt: on
    data_mask = np.where(data == no_data, np.NaN, data)
    driver = gdal.GetDriverByName(file_type)
    ds = driver.Create(filepath, pixel_x, pixel_y, num_band, dtype)

    if num_band == 1:
        band = ds.GetRasterBand(1)
        band.WriteArray(data_mask)
        band.SetNoDataValue(np.NaN)
        band.FlushCache()
    else:
        for i in range(0, num_band):
            band = ds.GetRasterBand(i + 1)
            band.WriteArray(data_mask[i, :, :])
            band.SetNoDataValue(np.NaN)
            band.FlushCache()

    if geomode == "proj":
        ds.SetGeoTransform(geotrans)
        ds.SetProjection(geoproj)
    elif geomode == "gcp":
        ds.SetGCPs(gcps, gcpsrc)

    ds.SetMetadata({'TIFFTAG_COPYRIGHT': 'KN_TOOL'}) 
    ds.SetMetadata({'TIFFTAG_XRESOLUTION': '1/1'})
    ds.SetMetadata({'TIFFTAG_YRESOLUTION': '1/1'})

    ds.FlushCache()
    return ds


def make_raster_from_array_and_prop(data, filepath, prop, dtype=gdal.GDT_Float32, no_data=9.9E33, file_type="GTiff"):

    if len(data.shape) == 2:
        num_band = 1
    else:
        num_band = data.shape[0]

    pixel_x = prop[0]
    pixel_y = prop[1]
    geotrans = prop[2]
    geoproj = prop[3]

    make_raster_from_array(
        data,
        filepath,
        pixel_x,
        pixel_y,
        num_band,
        dtype,
        no_data,
        file_type,
        "proj",
        geotrans=geotrans,
        geoproj=geoproj,
    )


def make_raster_with_gcps_from_array(data, filepath, dtype, no_data, prop, file_type="GTiff"):

    if len(data.shape) == 2:
        num_band = 1
    else:
        num_band = data.shape[0]

    pixel_x = prop[0]
    pixel_y = prop[1]
    gcpsrc = prop[2]
    gcps = prop[3]

    for gcp in gcps:
        lon = gcp.GCPX
        if lon < 0:
            gcp.GCPX = lon + 360.0

    make_raster_from_array(
        data, filepath, pixel_x, pixel_y, num_band, dtype, no_data, file_type, "gcp", gcps=gcps, gcpsrc=gcpsrc
    )


def make_geoproj(epsg):
    geoproj = osr.SpatialReference()
    geoproj.ImportFromEPSG(epsg)
    geoproj_wkt = geoproj.ExportToWkt()
    return geoproj_wkt


def make_south_nsidc_geoinfo(res):
    geotrans = ((-3950) * 1000, res * 1000, 0.0, (4350) * 1000, 0.0, -res * 1000)
    geoproj = make_geoproj(3412)
    return geotrans, geoproj


def make_north_nsidc_geoinfo(res):
    geotrans = ((-3850) * 1000, res * 1000, 0.0, (5850) * 1000, 0.0, -res * 1000)
    geoproj = make_geoproj(3411)
    return geotrans, geoproj


def open_generic_binary(infile, band, length, width, in_dtype=np.float32, byte_order="little"):
    with open(infile, mode="rb") as f:
        data = np.fromfile(f, dtype=in_dtype, sep="").reshape(band, length, width)
    if byte_order == "big":
        data = data.byteswap()
    return data


def make_generic_binary(data, infile, in_dtype=np.float32, byte_order="little"):
    if byte_order == "big":
        data = data.byteswap()
    data.astype(in_dtype).tofile(infile)


def open_hdf(infile, in_dtype, out_dtype, band, length, width):
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    for i, val in enumerate(ds.GetSubDatasets()):
        name = val  # [0].split(':')[4]
        ds[name] = gdal.Open(ds.GetSubDatasets()[i][0], gdal.GA_ReadOnly)
    return ds


def make_south_nsidc_raster_from_polygon(infile, outfile, width, length, band, file_type,
                                         in_dtype, byte_order, out_dtype, res, no_data):
    source_ds = ogr.Open(infile)
    source_layer = source_ds.GetLayer()

    geotrans, geoproj = make_south_nsidc_geoinfo(res)

    driver = gdal.GetDriverByName(file_type)

    ds = driver.Create(outfile, width, length, band, out_dtype[1])
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(geoproj)

    band = ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.FlushCache()

    gdal.RasterizeLayer(ds, [1], source_layer, options=["ATTRIBUTE=id", "ALL_TOUCHED=TRUE"])  # , burn_values=[1]
    #    gdal.RasterizeLayer(ds, [1], source_layer,burn_values=[1])
    ds.FlushCache()


def export_vector_from_geom(infile, geometry, dict=None, epsg=4326, srs=None):

    ext = os.path.splitext(infile)[::-1][0]

    if ext == ".shp":
        ftype = "ESRI Shapefile"
    elif ext == ".kml":
        ftype = "KML"
    else:
        print("The file type is not supported currently.")
        return

    driver = ogr.GetDriverByName(ftype)
    ds = driver.CreateDataSource(infile)

    if srs is None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
    # create one layer

    layer = ds.CreateLayer("layer", srs, geometry.GetGeometryType())
    if dict is None:
        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(idField)
    # Create the feature and set values
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(geometry)
    # feature.SetField("id", 1)
    layer.CreateFeature(feature)
    feature = None
    # Save and close DataSource
    ds = None


def export_vector_from_geomList(infile, geom_list, attr_dict=None, epsg=4326, srs=None):

    ext = os.path.splitext(infile)[::-1][0]

    if ext == ".shp":
        ftype = "ESRI Shapefile"
    elif ext == ".kml":
        ftype = "KML"
    else:
        print("The file type is not supported currently.")
        return

    driver = ogr.GetDriverByName(ftype)
    ds = driver.CreateDataSource(infile)
    if srs is None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)

    layer = ds.CreateLayer("layer", srs, geom_list[0].GetGeometryType())
    if attr_dict is None:
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(idField)
    else:
        keys = attr_dict.keys()
        [layer.CreateField(ogr.FieldDefn(key, OGRTypes[type(attr_dict[key][0])])) for key in keys]

    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)

    for i, geom in enumerate(geom_list):
        feature.SetGeometry(geom)
        if attr_dict is not None:
            [feature.SetField(key, attr_dict[key][i]) for key in keys]
        layer.CreateFeature(feature)
    feature = None
    # Save and close DataSource
    ds = None
