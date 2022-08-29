# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:21:32 2021

@author: Yo
"""
import numpy as np
from tkinter import sys
from scipy import signal
from scipy.ndimage import median_filter
from scipy.stats import linregress, zscore
import soundfile as sf
import pandas as pd

import time 

def cut_mono(x, fs):
    '''
    Parameters
    ----------
    x : array
        Input RIR signal.
    fs : int
        RIR's sampling frequency.

    Returns
    -------
    x_cortada : array
        Trimmed RIR signal.

    '''
    
    # Windowing starting from the maximum value
    
    x_max = np.where(abs(x) == np.max(abs(x)))[0][0]
    x_cortada = x[(x_max)+5:]
    
    # Trim the signal down to 10 seconds if the file is excessively long
    
    if len(x_cortada) / fs > 10: 
        x_cortada = x_cortada[:int(10 * fs)]    
       
    return x_cortada   

def cut_stereo(L, R, fs):
    '''
    Parameters
    ----------
    L : array
        Left channel RIR.
    R : array
        Right channel RIR.
    fs : int
        Sampling frequency.

    Returns
    -------
    L_cortada : array
        Trimmed left channel signal.
    R_cortada : TYPE
        Trimmed right channel signal.

    '''
    
    # Trim both signals
    
    L_cortada = cut_mono(L, fs)
    R_cortada = cut_mono(R, fs)
    
    # Get both signals to equal length
    
    if R_cortada.size > L_cortada.size:
        R_cortada[:L_cortada.size]    
    else:
        L_cortada[:R_cortada.size]
    
    return L_cortada, R_cortada

    
def filtroter(IR_L, IR_R, fs, ter):
    '''
    Filter a stereo impulse response according to the UNE-EN 61260 standard,
    using passband Butterworth filters.
    
    Parameters
    ----------
    IR_L : array
        Left channel input signal.
    IR_R : array
        Right channel input signal.
    fs : int
        Sampling frequency.
    ter : bool
        Filter by third octave (True) or octave band (False).

    Returns
    -------
    IR_L_filt : array
        Filtered left channel impulse response.
    IR_R_filt : array
        Filtered right channel impulse response.
    centrosHZ : array
        Filter's center frequencies.

    '''
    
    filtradaL = filtroter_mono(IR_L, fs, ter)
    filtradaR = filtroter_mono(IR_R, fs, ter)
    
    IR_L_filt = filtradaL[0]
    IR_R_filt = filtradaR[0]
    centrosHZ = filtradaL[1]
                    
    return  IR_L_filt, IR_R_filt, centrosHZ

def filtroter_mono(IR, fs, ter):
    '''
    Filter an impulse response according to the UNE-EN 61260 standard, using
    passband Butterworth filters.   
    
    Parameters
    ----------
    IR : array
        Input signal.
    fs : int
        Sampling frequency.
    ter : bool
        Filter by third octave (True) or octave band (False).

    Returns
    -------
    IR_filt : array
        Filtered impulse response.
    centrosHZ : array
        Filter's center frequencies.

    '''
    
    # Invert the impulse response and initialize variables
    
    W = np.flip(IR) 
    G = 10**(3/10)
    fil = []

    if ter:
        centrosHZ = np.array([25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 
                              250, 315, 400, 500, 630, 800, 1000, 1250, 1600,
                              2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000,
                              12500, 16000, 20000])
        fmin = G ** (-1/6)
        fmax = G ** (1/6)
    else:
        centrosHZ = np.array([31.5, 63, 125, 250, 500, 1000, 2000, 4000,
                              8000, 16000])
        fmin = G ** (-1/2)
        fmax = G ** (1/2)
            
    for j, fc in enumerate(centrosHZ):
        
        # Define the upper and lower limits of the frequency band

        sup = fmax*fc/(0.5*fs) 
        
        if sup >= 1:
            sup = 0.999999
            
        inf = fmin * fc / (0.5*fs) # Límite inferior

    # Apply the Nth order IIR Butterworth filter.
        
        sos = signal.butter(N=2, Wn=np.array([inf, sup]), 
                            btype='bandpass',output='sos')
        
        filt = signal.sosfilt(sos, W)
        fil.append(filt) 
        fil[j] = np.flip(fil[j])
    
    IR_filtrada = np.array(fil)
    
    # Cut the last 5% of the signal to minimize the border effect
    
    IR_filt = IR_filtrada[:int(len(fil[1])*0.95)]
    
    return IR_filt, centrosHZ

def ttyedtt(x, y, fs):
    '''
    Calculate Tt and early decay time (EDT) from a filtered IR.

    Parameters
    ----------
    x : array
        Array containing the filtered IR per frequency band in each row.
    y : array
        Array containing the .
    fs : int
        Sampling frequency.

    Returns
    -------
    Tt : float
        Transition time.
    EDTt : float
        Early decay time.

    '''
    
    EDTt = []
    Tt = []
    
    for i, ir in enumerate(x):
        
        # Remove the first 5 ms
        
        ir = ir[int(5e-3 * fs):]  
        
        # Find the index of the upper limit of the interval that contains 
        # 99% of the energy 
        
        index = np.where(np.cumsum(ir ** 2) <= 0.99 * np.sum(ir ** 2))[0][-1]
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

def mediana_movil(x, vent, fs, p):
    '''
    Filters the input signal x using a moving median filter with a window size
    of vent

    Parameters
    ----------
    x : array
        Array containing the input signal.
    vent : float
        Window size in milliseconds.
    fs : int
        Sampling frequency.
    p : int
        DESCRIPTION.

    Returns
    -------
    filt : array
        Filtered signal.

    '''

    v = int(vent*fs/1000)
    if v % 2 == 0:
        v +=1
    
    filt = median_filter(x,v)
    
    if p:
        filt = np.concatenate((filt, np.zeros(p)))

    with np.errstate(divide='ignore', invalid='ignore'):
        filt = 10*np.log10(filt / np.max(filt))
    return filt

def lundeby(x, fs):
    '''
    Finds the limits of integration of the Schroeder inverse integral
    according to Lundeby's method.

    Parameters
    ----------
    x : array
        Array containing the RIR signal.
    fs : int
        Sampling frequency.

    Returns
    -------
    punto : int
        Index of the upper limit of the Schroeder integral.
    C : float
        .

    '''

    N = x.size
    energia = x
    media = np.zeros(int(N/(fs*0.01)))
    eje_tiempo = np.zeros(int(N/(fs*0.01)))
    
    # Divide in sections and calculate the mean
    
    t = np.floor(N/(fs*0.01)).astype('int')
    v = np.floor(N/t).astype('int')
    
    for i in range(0, t):
        media[i] = np.mean(energia[i * v:(i + 1) * v])
        eje_tiempo[i] = np.ceil(v/2).astype('int') + (i*v)
        
    # Calculate noise level of the last 10% of the signal
    
    rms_dB = 10 * np.log10(np.sum(energia[round(0.9 * N):]) / (0.1 * N) / max(energia))
    mediadB = 10 * np.log10(media / max(energia))

    # Se busca la regresión lineal del intervalo de 0dB y la media mas proxima al ruido + 10dB.
    # Calculate linear regression between the 0 dB 
    
    try:
        r = int(max(np.argwhere(mediadB > rms_dB + 10)))
           
        if np.any(mediadB[0:r] < rms_dB+10):
            r = min(min(np.where(mediadB[0:r] < rms_dB + 10)))
        if np.all(r==0) or r<10:
            r=10
    except:
        r = 10

    # Least squares
        
    A = np.vstack([eje_tiempo[0:r], np.ones(len(eje_tiempo[0:r]))]).T
    m, c = np.linalg.lstsq(A, mediadB[0:r], rcond=-1)[0]
    cruce = int((rms_dB-c)/m)
    
    # Insufficient SNR
    
    if rms_dB > -20:
        
        punto = len(energia)
        C = None
        
    else:

        error = 1
        INTMAX = 50
        veces = 1
               
        while error > 0.0004 and veces <= INTMAX:
            
            # Calculates new time intervals for the mean with approximately
            # p steps for each 10 dB
            
            p = 10
            
            # Number of samples for the decay slope of 10 dB
            
            delta = int(abs(10/m)) 
            
            # Interval over which the mean is calculated
            
            v = np.floor(delta/p).astype('int') 
            t = int(np.floor(len(energia[:int(cruce-delta)])/v))
            
            if t < 2:
                t = 2
            elif np.all(t == 0):
                t = 2

            media = np.zeros(t)
            eje_tiempo = np.zeros(t)
            
            for i in range(0, t):
                media[i] = np.mean(energia[i*v:(i + 1) * v])
                eje_tiempo[i] = np.ceil(v / 2) + (i * v).astype('int')
                
            mediadB = 10 * np.log10(media / max(energia))
            A = np.vstack([eje_tiempo, np.ones(len(eje_tiempo))]).T
            m, c = np.linalg.lstsq(A, mediadB, rcond=-1)[0]

            # Nueva media de energia de ruido, comenzando desde desde el punto de la linea de 
            # decaimiento, 10 dB por debajo del punto de cruce
            
            # New noise average level, starting from the point of the 
            # decay curve, 10 dB below the intersection.
            
            noise = energia[int(abs(cruce + delta)):]
            
            if len(noise) < round(0.1 * len(energia)):
                noise = energia[round(0.9 * len(energia)):]
                
            rms_dB = 10 * np.log10(sum(noise)/ len(noise) / max(energia))

            # New intersection index
            
            error = abs(cruce - (rms_dB - c) / m) / cruce
            cruce = round((rms_dB - c) / m)
            veces += 1
                   
    # Output validation
            
    if cruce > N:
        punto = N
    else:
        punto = int(cruce)
        
    C = max(energia) * 10 ** (c / 10) * np.exp(m/10/np.log10(np.exp(1))*cruce) / (
        -m / 10 / np.log10(np.exp(1)))
        
    return punto, C

def RT_parameters(filtered_irs, fs):
    '''
    Calculate reverberation time parameters: Early Decay Time, T20 and T30.

    Parameters
    ----------
    filtered_irs : array
        Array containing the filtered IR per frequency band in each row.
    fs : int
        Frequency band.

    Returns
    -------
    EDT : float
    T20 : float
    T30 : float

    '''
    
    # Initialize variables and create time array
    
    results = {'EDT': [],
               'T20': [],
               'T30': []}
   
    t = np.arange(0, len(filtered_irs[0])/fs, 1/fs)
    
    for ir in filtered_irs:
    
        # Look for maximum values and their indexes
        
        i_max = np.where(ir == max(ir))[0][0]           
        y = ir[int(i_max):]
        y_max = max(y)
        
        # Get the indexes where the level of the signal is between the
        # defined limits
        
        i_edt = np.where((y <= y_max - 1) & (y > (y_max - 10)))
        i_20 = np.where((y <= y_max - 5) & (y > (y_max - 25)))    
        i_30 = np.where((y <= y_max - 5) & (y > (y_max - 35)))
        
        indexes = [i_edt, i_20, i_30]
        
        for i, key in zip(indexes, results):
            
            t_tr = np.vstack([t[i], np.ones(len(t[i]))]).T
            y_tr = y[i]
            
            # Calculate linear regression to extrapolate reverberation time
            # and append results to dictionary.
            
            m, c = np.linalg.lstsq(t_tr, y_tr, rcond=-1)[0]
            result = -60/m
            results[key].append(result)
        
    EDT = np.round(results['EDT'], 2)
    T20 = np.round(results['T20'], 2)
    T30 = np.round(results['T30'], 2)
    
    return EDT, T20, T30

def C_parameters(filtered_IRs, fs, lb):
    '''
    Calculates the clarity parameters C50 and C80 for a given RIR.

    Parameters
    ----------
    filtered_IRs : array
        Array containing the filtered RIR for each frequency band in each row.
    fs : int
        Sampling frequency.
    lb : array
        Array containing the limits of integration.

    Returns
    -------
    C50 : float
    C80 : float

    '''
    
    C50 = []
    C80 = []   
    
    for y, z in zip(filtered_IRs, lb):
        t50 = np.int64(0.05*fs)  # Index of signal value at 50 ms
        t80 = np.int64(0.08*fs)  # Index of signal value at 80 ms
        
        y50_num = y[0:t50]
        y50_den = y[t50:z]
        y80_num = y[0:t80]
        y80_den = y[t80:z]
    
        c50 = 10*np.log10(np.sum(y50_num) / np.sum(y50_den))
        c80 = 10*np.log10(np.sum(y80_num) / np.sum(y80_den))
        C50.append(c50)
        C80.append(c80)
        
    C50 = np.round(C50,2)
    C80 = np.round(C80,2)
    
    return C50, C80


def schroeder(x, p):
    '''
    Calculates the inverse Schroeder integral over an input RIR x.

    Parameters
    ----------
    x : array
        Array containing the RIR.
    p : int
        DESCRIPTION.

    Returns
    -------
    sch_db : array
        Inverse Schroeder integral curve in dB.

    '''
    
    sch = np.cumsum(x[::-1])[::-1]
    
    if p:
        sch = np.concatenate((sch, np.zeros(p)))

    with np.errstate(divide='ignore', invalid='ignore'):
        sch_db = 10.0 * np.log10(sch / np.max(sch))
    
    return sch_db

def E_norm(x):
    '''
    Calculates the energy time curve for a given input signal x.

    Parameters
    ----------
    x : array
        Input signal.

    Returns
    -------
    ETC : array
        Energy time curve.

    '''
    
    ETC = np.zeros(x.shape)  
    
    for i, y in enumerate(x):
        E = np.abs(signal.hilbert(y))**2
        ETC[i] = E/np.max(E)
        
    return ETC

def IACC_e(L, R, fs):
    '''
    Calculate IACCe according to the ISO 3382:2001 standard.

    Parameters
    ----------
    L : array
        Left channel input RIR.
    R : array
        Right channel input RIR.
    fs : int
        Sampling frequency.

    Returns
    -------
    IACCe : float
        Early interaural cross-correlation coefficient parameter.

    '''

    IACCe = []
    
    for ir_L, ir_R in zip(L, R):
        t80 = np.int64(0.08*fs)
        I = np.correlate(ir_L[0:t80], ir_R[0:t80], 'full')/(np.sqrt(np.sum(ir_L[0:t80]**2)*np.sum(ir_R[0:t80]**2)))
        iacce = np.max(np.abs(I))
        
        IACCe.append(iacce)
        
    IACCe = np.round(IACCe, 2)
        
    return IACCe

def ac_parameters_mono(W, ter, trunc, smooth, vent, fs):
    '''
    Calculates acoustical parameters from a single RIR.

    Parameters
    ----------
    W : array
        Input RIR signal.
    ter : bool
        Third octave (True) or octave (False) filter.
    trunc : bool
        Apply Lundeby (True) or not (False).
    smooth : int
        Schroeder (0), mmf (1) or no smoothing.
    vent : float
        Window size in ms.
    fs : int
        Sampling Frequency.

    Returns
    -------
    params : dict
        Acoustical parameters (Tt, EDTt, C50, C80, EDT, T20, T30).

    '''
    
    W = cut_mono(W, fs)
    
    # Filter RIR
    
    mono, cen = filtroter_mono(W, fs, ter)
    mono = E_norm(mono)
    
    # Obtain Energy Time Curve and smoothed curve

    etc = E_norm(np.array([W]))
    etc = etc[0]
    ETC = 10 * np.log10(etc + sys.float_info.epsilon)
        
    ir_tr = []
    ir_sm = []
    lb = []
      
    signals = np.append(mono, np.array([etc]), axis=0)
    
    for i, x in enumerate(signals):
        # Lundeby
        if trunc:
            ir_tr.append(x)
            punto_cruce = x.size
            lb.append(punto_cruce)
            p = 0
        else:
            punto_cruce, c = lundeby(x, fs)
            lb.append(punto_cruce)
            ir_tr.append(x[:punto_cruce])
            p = x.size-punto_cruce
        
        # Smoothing (Schroeder or MMF)
        if smooth == 1:
            mmf = mediana_movil(ir_tr[i], vent, fs, p)
            ir_sm.append(mmf)
            
        elif smooth == 0:
            sch = schroeder(ir_tr[i], p)
            ir_sm.append(sch)

    ir_tr.pop(-1)
    smooth_ir = ir_sm.pop(-1)
    
    Tt, EDTt = ttyedtt(ir_tr, ir_sm, fs)
    EDT, T20, T30 = RT_parameters(ir_sm, fs)
    C50, C80 = C_parameters(mono, fs, lb)
    
    params = {'Tt': Tt,
              'EDTt': EDTt,
              'C50': C50,
              'C80': C80,
              'EDT': EDT,
              'T20': T20,
              'T30': T30,
              'ETC': ETC,
              'smooth': smooth_ir
              }
    
    return params

def ac_parameters_stereo(L, R, ter, trunc, smooth, vent, fs):
    '''
    Calculates acoustical parameters from a stereo RIR.

    Parameters
    ----------
    L : array
        Input left RIR signal.
    R : array
        Input right RIR signal.
    ter : bool
        Third octave (True) or octave (False) filter.
    trunc : bool
        Apply Lundeby (True) or not (False).
    smooth : int
        Schroeder (0), mmf (1) or no smoothing.
    vent : float
        Window size in ms.
    fs : int
        Sampling Frequency.

    Returns
    -------
    params : dict
        Acoustical parameters (Tt, EDTt, C50, C80, EDT, T20, T30, IACCe).

    '''
    
    IR_L, IR_R = cut_stereo(L, R, fs)
    L_fil, R_fil, _ = filtroter(IR_L, IR_R, fs, ter)
    
    # Add IACCe parameter to results and append curves of both channels
    
    IACCe = IACC_e(L_fil, R_fil, fs)
    
    paramsL = ac_parameters_mono(L, ter, trunc, smooth, vent, fs)
    paramsR = ac_parameters_mono(R, ter, trunc, smooth, vent, fs)
    
    params = paramsL.copy()
    
    params['IACCe'] = IACCe
    params['ETC'] = [paramsL['ETC'], paramsR['ETC']]
    params['smooth'] = [paramsL['smooth'], paramsR['smooth']]
    
    return params

def create_table(params, ter):
    '''
    Creates a pandas dataframe containing the acoustical parameters.

    Parameters
    ----------
    params : dict
        Acoustical parameters.
    ter : bool
        Third octave (True) or octave (False) filter.

    Returns
    -------
    df : dataframe
        Pandas dataframe containing the acoustical parameters.

    '''
    
    data = params.copy()
    
    if ter:
        freqs = ['25', '31.5', '40', '50', '63', '80', '100', '125', '160',
                 '200', '250', '315', '400', '500', '630', '800', '1k','1.3k',
                 '1.6k', '2k', '2.5k', '3.2k', '4k', '5k', '6.3k', '8k', '10k', 
                 '12.5k', '16k', '20k']
    else:
        freqs = ['31.5', '63', '125', '250', '500', '1000', 
                 '2000', '4000', '8000', '16000']
    
    # Create DataFrame
    
    del data['ETC']
    del data['smooth']
    
    df = pd.DataFrame.from_dict(data, 
                                orient='index',
                                columns=freqs)
   
    return df

#%% PRUEBAS
    
# STEREO    
    
# pathL = 'IR S1 Posición 1 Mic Kemar L.wav'
# pathR = 'IR S1 Posición 1 Mic Kemar R.wav'
# dataL, fs = sf.read(pathL)
# dataR, fs = sf.read(pathR)

# L = np.transpose(dataL)
# R = np.transpose(dataR)

# divi = 0 #Filtro de octava (0) o tercio de octava (1)
# trunc = 0 #Lundeby (0) o None (1)
# smooth = 0 #Schroeder (0) o mmf (1)
# vent = 20 #Tamaño de ventana mmf

# params = ac_parameters_stereo(L, R, divi, trunc, smooth, vent, fs)
# data = create_table(params, divi)

# MONO

# start = time.time()

# W, fs = sf.read('IR S1 Posición 1 Mic Earthworks 1 (1).wav')

# divi = 0 #Filtro de octava (0) o tercio de octava (1)
# trunc = 0 #Lundeby (0) o None (1)
# smooth = 0 #Schroeder (0) o mmf (1)
# vent = 20 #Tamaño de ventana mmf
# params = ac_parameters_mono(W, divi, trunc, smooth, vent, fs)
# data = create_table(params, divi)

# end = time.time()

# print(end-start)
        
        
    
