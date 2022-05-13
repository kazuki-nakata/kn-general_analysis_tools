import arcpy

inlayer=r'inpath'
outlayer=r'outpath'
arcpy.RasterToPolygon_conversion(inlayer,outlayer,"SIMPLIFY","VALUE")