from ctypes import sizeof
from turtle import shape, shapesize
import scipy.signal as sc
import scipy.ndimage as scn
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
import acoustics.room as ac
import wavfile

# data,fs= sf.read(r"IRpromedio.wav",always_2d=True)

def tercio(data,fs):
    indice = np.arange(-16,13,1) # int numbers array [-16,-15,...,12,13]
    fr = 1000 # 
    b = 3
    fm_v = []
    f1_v = []
    f2_v = []

    signals = []

    for i in range(0,len(indice)):
        fm = fr*2**(indice[i]/b)
        f1 = fm*(2**(-1/(2*b)))/(0.5*fs)
        f2 = fm*(2**(1/(2*b)))/(0.5*fs)
        fm_v.append(fm)
        f1_v.append(f1)
        f2_v.append(f2)

    for i in range(0,len(indice)):
        sos = sc.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',output='sos')
        sig_filtrada = sc.sosfilt(sos, data)
        signals.append(sig_filtrada)
        
    return signals

def t60_impulse(file_name, bands, rt='t30'):  # pylint: disable=too-many-locals
    """
    Reverberation time from a WAV impulse response.

    :param file_name: name of the WAV file containing the impulse response.
    :param bands: Octave or third bands as NumPy array.
    :param rt: Reverberation time estimator. It accepts `'t30'`, `'t20'`, `'t10'` and `'edt'`.
    :returns: Reverberation time :math:`T_{60}`

    """
    raw_signal,fs = sf.read(file_name)
    band_type = ac._check_band_type(bands)

    if band_type == 'octave':
        low = ac.octave_low(bands[0], bands[-1])
        high = ac.octave_high(bands[0], bands[-1])
    elif band_type == 'third':
        low = ac.third_low(bands[0], bands[-1])
        high = ac.third_high(bands[0], bands[-1])

    rt = rt.lower()
    if rt == 't30':
        init = -5.0
        end = -35.0
        factor = 2.0
    elif rt == 't20':
        init = -5.0
        end = -25.0
        factor = 3.0
    elif rt == 't10':
        init = -5.0
        end = -15.0
        factor = 6.0
    elif rt == 'edt':
        init = 0.0
        end = -10.0
        factor = 6.0

    t60 = np.zeros(bands.size)

    for band in range(bands.size):
        # Filtering signal
        filtered_signal = ac.bandpass(raw_signal, low[band], high[band], fs, order=8)
        abs_signal = np.abs(filtered_signal) / np.max(np.abs(filtered_signal))

        # Schroeder integration
        sch = np.cumsum(abs_signal[::-1]**2)[::-1]
        sch_db = 10.0 * np.log10(sch / np.max(sch))

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
        t60[band] = factor * (db_regress_end - db_regress_init)
    return t60

bandas = [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,10000, 12500, 16000]
bandas= np.array(bandas)

t20=t60_impulse(r"IRpromedio.wav",bandas,rt='t20')

print(t20)

# signals=tercio(data.T,fs)
# signals=np.squeeze(signals)
# signals=abs(signals).flatten()
# signals = abs(signals)
# print("filtro por tercio",signals)
# sf.write("tercio.wav",signals,fs)

# maksimo=np.max(signals)

# print(maksimo)
##Filtro de Media movil
# test=scn.median_filter(signals,size=3)
# print("TESTMA:   ",len(test))
# print(type(test))
# def E_norm(x):
#     '''
#     Calculates the energy time curve for a given input signal x.
#     Parameters
#     ----------
#     x : array
#         Input signal.
#     Returns
#     -------
#     ETC : array
#         Energy time curve.
#     '''
#     dict(enumerate(x.flatten(), 1))
#     ETC = np.zeros(x.shape)  
    
#     for i, y in enumerate(x):
#         E = np.abs(sc.hilbert(y))**2
#         ETC[i] = E/np.max(E)
        
#     return ETC

# ETC=E_norm(signals)

# normalizeta=20*np.log10(signals/(np.max(maksimo)))
# print(normalizeta)


# plt.plot(ETC)
# plt.show()

bandas= [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,10000, 12500, 16000]
bandas=np.array(bandas)

# EDT,T20,T30 =RT_parameters(test,fs)

# def t60_impulse(file_name, bands, rt='t30'):  # pylint: disable=too-many-locals
#     """
#     Reverberation time from a WAV impulse response.

#     :param file_name: name of the WAV file containing the impulse response.
#     :param bands: Octave or third bands as NumPy array.
#     :param rt: Reverberation time estimator. It accepts `'t30'`, `'t20'`, `'t10'` and `'edt'`.
#     :returns: Reverberation time :math:`T_{60}`

#     """
#     fs, raw_signal = wavfile.read(file_name)
#     band_type = ac._check_band_type(bands)

#     if band_type == 'octave':
#         low = ac.octave_low(bands[0], bands[-1])
#         high = ac.octave_high(bands[0], bands[-1])
#     elif band_type == 'third':
#         low = ac.third_low(bands[0], bands[-1])
#         high = ac.third_high(bands[0], bands[-1])

#     rt = rt.lower()
#     if rt == 't30':
#         init = -5.0
#         end = -35.0
#         factor = 2.0
#     elif rt == 't20':
#         init = -5.0
#         end = -25.0
#         factor = 3.0
#     elif rt == 't10':
#         init = -5.0
#         end = -15.0
#         factor = 6.0
#     elif rt == 'edt':
#         init = 0.0
#         end = -10.0
#         factor = 6.0

#     t60 = np.zeros(bands.size)

#     for band in range(bands.size):
#         # Filtering signal
#         filtered_signal = bandpass(raw_signal, low[band], high[band], fs, order=8)
#         abs_signal = np.abs(filtered_signal) / np.max(np.abs(filtered_signal))

#         # Schroeder integration
#         sch = np.cumsum(abs_signal[::-1]**2)[::-1]
#         sch_db = 10.0 * np.log10(sch / np.max(sch))

#         # Linear regression
#         sch_init = sch_db[np.abs(sch_db - init).argmin()]
#         sch_end = sch_db[np.abs(sch_db - end).argmin()]
#         init_sample = np.where(sch_db == sch_init)[0][0]
#         end_sample = np.where(sch_db == sch_end)[0][0]
#         x = np.arange(init_sample, end_sample + 1) / fs
#         y = sch_db[init_sample:end_sample + 1]
#         slope, intercept = stats.linregress(x, y)[0:2]

#         # Reverberation time (T30, T20, T10 or EDT)
#         db_regress_init = (init - intercept) / slope
#         db_regress_end = (end - intercept) / slope
#         t60[band] = factor * (db_regress_end - db_regress_init)
#     return t60


# test=test/np.max(test)

# Tparam=ac.t60_impulse(r"IRpromedio.wav",bandas,rt='t20')
# Tparam=t60_impulse(np.array(test),bandas,fs,rt='t20')

# Tparam=20*np.log10(Tparam)
# print(Tparam)
# plt.stem(Tparam)
# plt.show()
# signals=np.array(signals)
# fsignals=media_misley(signals,3)
# print("MEDIA MOVIL:   ",fsignals)

# sf.write("media.wav",signals,fs)
# print(arsignals.shape)

## Time plots
# plt.figure()
# plt.subplot(211)
# plt.plot(test)
# plt.title('Test')
# plt.subplot(212)
# plt.plot(signals)
# plt.title('Media misley')
# plt.show()

# signals=np.fft.fft(signals)
# fsignals=np.fft.fft(test)
# xfrecs=np.int32(len(signals))
# frecs=np.fft.fftfreq(xfrecs,1/fs)

# half = len(signals)//2
# signals=signals[:half]

#dB
# signals=20*np.log10(signals)
# fsignals=20*np.log10(fsignals)
#Full scale
# signals=20*np.log10(signals/np.max(signals))
# fsignals=20*np.log10(fsignals/np.max(fsignals))

# plt.figure()
# plt.subplot(211)
# plt.plot(frecs[:len(frecs)//2],signals[:len(frecs)//2])
# plt.title('filtro tercio')
# plt.subplot(212)
# plt.plot(frecs[:len(frecs)//2],fsignals[:len(frecs)//2])
# plt.title('Media movil')
# plt.show()

# signals2=tercio(data2,fs)
# arsignals2=np.array(signals2)
# fsignals=media_misley(arsignals2,10)
# print(arsignals2)
# print(arsignals2.shape)


# print(arsignals.all()==arsignals2.all())

# plt.subplot(211)
# plt.plot(arsignals[0])
# plt.subplot(212)
# plt.plot(arsignals2[0])
# plt.show()
# for signali in arsignals2:
#     plt.plot(arsignals2[signali])
    
# plt.show()
# print(arsignals)
# for signal in arsignals:
#     plt.plot(signal)
#     plt.show()
    


# datastereo=wavfile.read("RuidoRosa.wav",fmt='float')
# datastereo=np.array(datastereo)


# print(datastereo)