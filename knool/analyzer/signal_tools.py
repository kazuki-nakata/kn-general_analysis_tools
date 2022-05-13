import numpy as np
from scipy import signal
from scipy.fft import fft, fftfreq #not recommend fftpack
from scipy.signal import hilbert
import pandas as pd

#https://www.logical-arts.jp/archives/112
def fft_overlap(data, samplerate, Fs, overlap):
    Ts = len(data) / samplerate         #全データ長
    Fc = Fs / samplerate                #フレーム周期
    x_ol = Fs * (1 - (overlap/100))     #オーバーラップ時のフレームずらし幅
    N_ave = int((Ts - (Fc * (overlap/100))) / (Fc * (1-(overlap/100)))) #抽出するフレーム数（平均化に使うデータ個数）
 
    array = []      #抽出したデータを入れる空配列の定義
 
    #forループでデータを抽出
    for i in range(N_ave):
        ps = int(x_ol * i)              #切り出し位置をループ毎に更新
        array.append(data[ps:ps+Fs:1])  #切り出し位置psからフレームサイズ分抽出して配列に追加
    return array, N_ave                 #オーバーラップ抽出されたデータ配列とデータ個数を戻り値にする


def fft(data_array,samplerate, Fs, N_ave, acf):
    fft_array = []
    for i in range(N_ave):
        fft_array.append(acf*np.abs(fft(data_array[i])/(Fs/2))) #FFTをして配列に追加、窓関数補正値をかけ、(Fs/2)の正規化を実施。
    fft_axis = np.linspace(0, samplerate, Fs)   #周波数軸を作成
    fft_array = np.array(fft_array)             #型をndarrayに変換
    fft_mean = np.sqrt(np.mean(fft_array ** 2, axis=0))       #全てのFFT波形のパワー平均を計算してから振幅値とする
    return fft_array, fft_mean, fft_axis

def window_function_list(csort=None):
    #https://org-technology.com/posts/scipy-window-function.html
    col_name=['name','side_lev(dB)','main_width_3dB','main_width_6dB','Scallop Loss','NENBW']
    df=pd.read_csv(r'./window_function.csv',sep='\t',names=col_name)
    
    if csort != None:
        df=df.sort_values(csort)
    print(df)

def calc_shift_from_cc_1D(input1,input2):
    csd = np.fft.fft(input1) * np.fft.fft(input2).conj()
    output = np.fft.fftshift(np.fft.ifft(csd))
    return output

def calc_shift_from_cc_2D(input1,input2):
    csd = np.fft.fft2(input1) * np.fft.fft2(input2).conj()
    output = np.fft.fftshift(np.fft.ifft2(csd))
    return output

def hilbert_transform(inArray,win=1):
    outArray = hilbert(inArray*win)/win
    return outArray

def bandpass_filter(inArray, numtaps, cutoff, fs,pass_zero, window='hamming'): #just wrap
    fir_filter = signal.firwin(numtaps=numtaps,
                               cutoff=cutoff,
                               fs=fs,
                               pass_zero=pass_zero,
                               window=window)
    outArray = signal.lfilter(fir_filter,1,inArray)
    return outArray

