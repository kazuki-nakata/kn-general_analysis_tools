import sys
import numpy as np
from osgeo import osr

class Projection:
    def __init__(self, source_ref, target_ref,coord_array,index_array,res):
        #source_ref and target_ref: ref=osr.SpatialReference() -> ref.ImportFromEPSG(epsg)
        self.target_ref=target_ref
        self.res=res
        self.build(source_ref,target_ref,coord_array,index_array,res)
        
    def build(self,source_ref,target_ref,coord_array,index_array,res):
        #coord_array : For latlon, array(n,2). (n,0)->lat (n,1)->lon
        #output trans_array: array(n,2). (n,0) -> x axis (horizontal)
        coord_transform=osr.CoordinateTransformation(source_ref, target_ref)                        
        self.trans_array=np.array(coord_transform.TransformPoints(coord_array))[:,0:2]
        print("transform")
        print(coord_array[0,0:2],coord_array.shape)
        print(self.trans_array[0,0:2],self.trans_array.shape)
        
        self.trans_array=np.int32(self.trans_array/res)
        imin=np.min(self.trans_array[:,0])
        imax=np.max(self.trans_array[:,0])
        jmin=np.min(self.trans_array[:,1])
        jmax=np.max(self.trans_array[:,1])        
        self.nx=imax-imin+1
        self.ny=jmax-jmin+1
        self.offset_i=imin
        self.offset_j=jmin
        #------i -> horizon
        self.trans_array[:,0]=self.trans_array[:,0]-self.offset_i
        self.trans_array[:,1]=self.trans_array[:,1]-self.offset_j
                
        print(self.nx,self.ny,self.offset_i,self.offset_j)
        self.map_field_i=np.full([self.nx,self.ny],-1)
        self.map_field_j=np.full([self.nx,self.ny],-1)
        print(self.map_field_i.shape,index_array.shape,self.trans_array.shape)
       
        #should change array(n,2) to array(2,n)
        trans_array_t=self.trans_array.T
        
        print(trans_array_t)
        self.map_field_i[trans_array_t.tolist()]=index_array[0]
        self.map_field_j[trans_array_t.tolist()]=index_array[1]

        
    def search_index(self,source_ref,coord_array,index_array):
        #coord_array : For latlon, array(n,2). (n,0)->lat (n,1)->lon
        #output coord_array: array(n,2). (n,0) -> x axis (horizontal)
        coord_transform=osr.CoordinateTransformation(source_ref, self.target_ref)          
        trans_array=np.array(coord_transform.TransformPoints(coord_array))[:,0:2]
        trans_array=np.int32(trans_array/self.res)
        trans_array[:,0]=trans_array[:,0]-self.offset_i
        trans_array[:,1]=trans_array[:,1]-self.offset_j
        
        print(index_array.shape,trans_array.shape)
        trans_array2=trans_array[(trans_array[:,0]>0) & (trans_array[:,0]<self.nx) & (trans_array[:,1]>0) & (trans_array[:,1]<self.ny)].T
        index_array2=index_array[(trans_array[:,0]>0) & (trans_array[:,0]<self.nx) & (trans_array[:,1]>0) & (trans_array[:,1]<self.ny)].T
        print(index_array2.shape,trans_array2.shape)
        out_ij=np.array([self.map_field_i[trans_array2.tolist()],self.map_field_j[trans_array2.tolist()]])
        bool_array=(out_ij[0,:]!=-1)
        out_ij=out_ij[:,bool_array]
        index_array2=index_array2[:,bool_array]
        
        return out_ij,index_array2
    