import numpy as np
import scipy
from scipy import signal
from scipy.signal import hilbert
from scipy import hamming


def reorder_3d_to_2d(array,index): #indexはリスト。2番目と3番目で指定されたインデックスが結合される。
    shape=array.shape
    z,y,x=index
    zmax,ymax,xmax=shape[z],shape[y],shape[x]
    array2=array.transpose(z,y,x)
    array3=array2.reshape([zmax,ymax*xmax],order="C")
    return array3

def reorder_2d_to_3d(array,index,shape): #indexはリスト。shapeで２，３番目で指定された形で分解される。
    zmax,ymax,xmax=shape[0],shape[1],shape[2]
    array2=array.reshape([zmax,ymax,xmax],order="C")
    z,y,x=index
    array3=array2.transpose(z,y,x)
    return array3

def calc_PCA(array,mode="numpy",norm=True,limit=10):
    if mode=="numpy":
        cov_matrix=np.cov(array,rowvar=False) #rowvarは軸の切り替え
        eig_val, eig_vec = np.linalg.eigh(cov_matrix) 
        #eigは全固有値問題、eighは実対称もしくはエルミート行列のみ.出力はeig_vec[parameter,mode]
        #複素データに対し、eighで計算された固有ベクトルの最初の要素は必ず0orpi位相になる。（虚部が０）
    print(eig_val.shape,eig_vec.shape)
    eig_val_sorted=np.sort(eig_val)[::-1]
    eig_val_sortedi=np.argsort(eig_val)[::-1]
    eig_val_sum = np.sum(eig_val)
    for i in range(limit):
        cont = (eig_val_sorted[i] / eig_val_sum )* 100
        print("第",i+1,"主成分寄与率：",cont,"インデックス:",eig_val_sortedi[i])
    comCont = np.sum(eig_val_sorted[0:limit])/ np.sum(eig_val_sorted) * 100
    print("第",limit,"主成分までの累積寄与率:",comCont)
    return eig_val,eig_vec,eig_val_sortedi

def calc_PCA_score(array,eig_vec,index):
    score=np.dot(array, np.conj(eig_vec[:,index]))
    print(score.shape)
    return score

def reconst_data_from_score(array,eig_vec,mindex):
    score=calc_PCA_score(array,eig_vec,mindex)
    print(eig_vec[:,mindex].shape,np.conj(score).shape)
    result=np.dot(score,eig_vec[:,mindex].T)
    return result

def remove_noise_for_xyt(array,n_mode):
    tmax,ymax,xmax=array.shape
    reorder=reorder_3d_to_2d(array,[0,1,2])
    eig_val,eig_vec,sindex=calc_PCA(reorder)
    index=sindex[0:n_mode]
    reconst=reconst_data_from_score(reorder,eig_vec,index)
    result=reorder_2d_to_3d(reconst,[0,1,2],array.shape)
    return result

def calc_EOF_from_xyt(array,mmin,mmax):
    tmax,ymax,xmax=array.shape
    reorder=reorder_3d_to_2d(array,[0,1,2])
    eig_val,eig_vec,sindex=calc_PCA(reorder)
    eig_val_sorted=np.sort(eig_val)[::-1]
    mindex=sindex[mmin-1:mmax]
    temporal=calc_PCA_score(reorder,eig_vec,mindex)
    spatial=eig_vec[:,mindex].reshape([ymax,xmax,mmax-mmin+1]).transpose(2,0,1)
    return spatial,temporal,eig_val_sorted

