#%%
import numpy as np
import sys
import os
from ..helpers.misc import deco_print_time
from ..helpers.geometry import map_2D
import cython

#%%

@deco_print_time('start','end')
def calc_IST_with_split_window(bt1_array,bt2_array,senz_array,param): #bt1->low, bt2->high, senz->sensor zenith
    
    get_IST = lambda bt1, bt2, senz, p: p['a'] + p['b']*bt1 + p['c']*(bt1-bt2) + p['d']*(bt1-bt2)*(1./np.cos(np.deg2rad(senz))-1.)    
    outArray=np.where(bt1_array<240, get_IST(bt1_array,bt2_array,senz_array,param['240K']),
             np.where((bt1_array>=240)&(bt1_array<260), get_IST(bt1_array,bt2_array,senz_array,param['240-260K']),
             get_IST(bt1_array,bt2_array,senz_array,param['260K'])))
    
    return outArray

