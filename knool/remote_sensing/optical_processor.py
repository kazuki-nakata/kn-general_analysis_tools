import os, time, gdal
import numpy as np

def RGB_for_byte_stack(red_path, green_path, blue_path, prop, max, min, no_data):
    red = gdal.Open(red_path, gdal.GA_ReadOnly)
    green = gdal.Open(green_path, gdal.GA_ReadOnly)
    blue = gdal.Open(blue_path, gdal.GA_ReadOnly)

    red_arr = red.ReadAsArray().astype(np.float16)
    red_arr[red_arr == no_data] = np.nan
    green_arr = green.ReadAsArray().astype(np.float16)
    green_arr[green_arr == no_data] = np.nan
    blue_arr = blue.ReadAsArray().astype(np.float16)
    blue_arr[blue_arr == no_data] = np.nan

    red_arr[red_arr < min] = min
    red_arr[red_arr > max] = max 
    green_arr[green_arr < min] = min
    green_arr[green_arr > max] = max
    blue_arr[blue_arr < min] = min
    blue_arr[blue_arr > max] = max

    red_arr_nor = red_arr * (255 / max)
    green_arr_nor = green_arr * (255 / max)
    blue_arr_nor = blue_arr * (255 / max)

    stack = np.empty((3, prop[1], prop[0]), dtype = np.float16)
    stack[0,:,:] = red_arr_nor
    stack[1,:,:] = green_arr_nor
    stack[2,:,:] = blue_arr_nor

    del red, green, blue	

    return stack