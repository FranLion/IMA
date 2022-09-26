from PyQt5.QtWidgets import *
import scipy.signal as sc
import scipy.stats as st
import scipy.ndimage as scn
import acoustics.room as ac
import soundfile as sf
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from librosa import amplitude_to_db
import numpy as np
import sys
from filterbank import Filterbank


def filtrado(rir, fs, ter):
    params = {'fs' : fs,
              'bands' :[125, 250, 500, 1000, 2000, 4000, 8000, 16000],
              'bandsize' : 1,
              'order' : 4,
              'f_length': 16384,
              'power' : True}
    if ter:
        params['bands'] = [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000]
        params['bandsize'] = 3
    filterbank = Filterbank(**params)
    bands, centros = filterbank.apply(rir)  
    return bands, centros


def denoise(signal, fs, ventana = 0.05):
    eps = np.finfo(float).eps
    last_10 = int(signal.size * 0.9)
    N = int(ventana * fs) 
    env = sc.fftconvolve(abs(signal), np.ones(N)/N, mode='same')
    min_value = np.min(np.nan_to_num(env))
    noise_floor = np.nan_to_num(env[last_10:], nan=min_value).mean()
    noise_floor_idx = np.where(env<noise_floor)[0][0]
    return signal[:noise_floor_idx]


def rt_descriptors(signal, signal_raw, fs):
    # signal ==> filtrada, suavizada
    
    # rts -> param : init, end, factor
    rts = {'t30' : [-5.0, -35.0, 2.0],
           't20' : [-5.0, -25.0, 3.0],
           't10' : [-5.0, -15.0, 6.0],
           'edt' : [0.0, -10.0, 6.0]}
    params = {}
    
    for rt in rts:
        init, end, factor = rts[rt]
        sch_db = signal
        # Linear regression
        sch_init = sch_db[np.abs(sch_db - init).argmin()]
        sch_end = sch_db[np.abs(sch_db - end).argmin()]
        init_sample = np.where(sch_db == sch_init)[0][0]
        end_sample = np.where(sch_db == sch_end)[0][0]
        x = np.arange(init_sample, end_sample + 1) / fs
        y = sch_db[init_sample:end_sample + 1]
        slope, intercept = ac.stats.linregress(x, y)[0:2]
        # Reverberation time (T30, T20, T10 or EDT)
        db_regress_init = (init - intercept) / slope
        db_regress_end = (end - intercept) / slope
        param = factor * (db_regress_end - db_regress_init)
        params[rt] = param
        
    # clarity -> param : c50/c80
    clarities = {'C50' : 50, 'C80' : 80}
    for clarity in clarities:
        h2 = signal_raw**2.0
        time = clarities[clarity]
        t = int((time / 1000.0) * fs + 1) #Así venía el original
        #t = int(time / 1000.0) * fs 
        c = 10.0 * np.log10((np.sum(h2[:t]) / np.sum(h2[t:])) + sys.float_info.epsilon)
        params[clarity] = c
    return params


# def filtrado(IR, fs, ter):
#     '''
#     Filter an impulse response according to the UNE-EN 61260 standard, using
#     passband Butterworth filters.   
    
#     Parameters
#     ----------
#     IR : array
#         Input signal.
#     fs : int
#         Sampling frequency.
#     ter : bool
#         Filter by third octave (True) or octave band (False).
#     Returns
#     -------
#     IR_filt : array
#         Filtered impulse response.
#     centrosHZ : array
#         Filter's center frequencies.
#     '''
#     # Invert the impulse response and initialize variables
#     W = np.flip(IR)
#     # W= IR
#     G = 10**(3/10)
#     fil = []
#     if ter:
#         centrosHZ = np.array([25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000])
#         fmin = G ** (-1/6)
#         fmax = G ** (1/6)
#     else:
#         centrosHZ = np.array([125, 250, 500, 1000, 2000, 4000, 8000, 16000])
#         fmin = G ** (-1/2)
#         fmax = G ** (1/2) 
#     for j, fc in enumerate(centrosHZ):
#         # Define the upper and lower limits of the frequency band
#         sup = fmax*fc/(0.5*fs)  
#         if sup >= 1:
#             sup = 0.999999    
#         inf = fmin * fc / (0.5*fs) # Límite inferior
#     # Apply the Nth order IIR Butterworth filter. 
#         sos = sc.butter(N=4, Wn=np.array([inf, sup]), btype='bandpass',output='sos')
#         filt = sc.sosfilt(sos, W)
#         fil.append(filt) 
#         fil[j] = np.flip(fil[j])
#     IR_filtrada = np.array(fil)
#     # Cut the last 5% of the signal to minimize the border effect
#     IR_filt = IR_filtrada[:int(len(fil[1])*0.95)]
#     return np.squeeze(IR_filt), centrosHZ

def ttyedtt(x, y, fs):
    # x : filtrada suavizada y en dB
    
    EDTt, Tt = [], []
    for i, ir in enumerate(x):
        # Remove the first 5 ms
        ir = ir[int(5e-3 * fs):]  
        # Find the index of the upper limit of the interval that contains 
        # 99% of the energy 
        index = np.where(np.cumsum(ir * 2) <= 0.99 * np.sum(ir * 2))[0][-1]
        t_t = index/fs
        Tt.append(t_t)
        # Filter the impulse with the moving median filter in order to calculate
        # the parameters
        ir2 = y[i]
        ir2 = ir2[:index]  
        t_Tt = np.arange(0, index/fs, 1/fs)
        
        if len(t_Tt) > index:
            t_Tt = t_Tt[:index]
        
        # Calculate minimum squares
            
        A = np.vstack([t_Tt, np.ones(len(t_Tt))]).T
        m, c = np.linalg.lstsq(A, ir2, rcond=-1)[0]
        
        edt_t = -60/m
        EDTt.append(edt_t)
        
    EDTt = np.round(EDTt, 2)
    Tt = np.round(Tt, 2)
    
    return Tt, EDTt

