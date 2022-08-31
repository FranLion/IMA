import scipy as sc
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
""" 
This function filters the input signal in octave and 1/3 octave bands.
Reference:
    fr: reference frequency
    fs: sampling frequency
    sig: input signal
"""

# x = np.arange(10000)
# f = 1000
# sig = np.sin(2 * np.pi * f * x / fs)

# sig = np.random.normal(0,1,100) 

sig, fs = sf.read('tone.wav')


indice = np.arange(-16,13,1) # int numbers array [-16,-15,...,12,13]
fs = 44100 # 
fr = 1000 # 
b = 3
fm_v = []
f1_v = []
f2_v = []
h_v = []
# leq = []
# rms_banda = []

# sig = sig + 2
# sig_cal , fs = sf.read('Calibracion ruido de fondo.wav')
# sig , fs = sf.read('Ruido de fondo.wav')
# rms= np.sqrt((1/len(sig_cal))*np.sum(sig_cal**2))
# sig = sig/rms # Señal calibrada
xticks = [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,10000, 12500, 16000]
xlabels = ['25' ,'31.5', '40','50','63','80','100','125','160','200','250','315','400','500','630','800','1k','1.25k','1.6k','2k','2.5k','3.15k','4k','5k','6.3k','8k','10k','12.5k','16k']

signals = []

for i in range(0,len(indice)):
    fm = fr*2**(indice[i]/b) 
    f1 = fm*(2**(-1/(2*b)))/(0.5*fs)
    f2 = fm*(2**(1/(2*b)))/(0.5*fs)
    fm_v.append(fm)
    f1_v.append(f1)
    f2_v.append(f2)

for i in range(0,len(indice)):
    sos = sc.signal.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',output='sos')
    a,b = sc.signal.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',analog=True,output='ba')
    w, h = sc.signal.freqs(a,b,worN=20000)
    # h_v.append(h)
    # plt.semilogx(w*fs/2, 20 * np.log10(abs(h)))
    sig_filtrada = sc.signal.sosfilt(sos, sig)
    signals.append(sig_filtrada) 
    # rms_banda=np.sqrt((1/len(sig_filtrada))*np.sum(sig_filtrada**2))
    # leq_señal=10*np.log10((rms_banda/20e-6)**2)
    # leq.append(leq_señal)

# plt.xticks(xticks, xlabels, rotation=70)
# plt.xlim(20,20000)
# plt.ylim(-20, 2.5)
# plt.xlabel(r'$Frecuencia\ [Hz]$', fontsize=12)
# plt.ylabel(r'$Magnitud\ [dB]$', fontsize=13)
# plt.grid()
# plt.tight_layout()
# plt.savefig('rta_filtros.png')

n=17
plt.plot(sig)
plt.plot(signals[n])

sf.write('filteredTone.wav', signals[n], fs)

# freqs = np.arange(0,29,1)

# plt.figure()
# plt.bar(range(len(fm_v)),leq)
# plt.xticks(freqs,xlabels,rotation=70)
# plt.ylabel(r'$L_{eq}\ [dB_{SPL}]$', fontsize=13)
# plt.xlabel(r'$Frecuencia\ [Hz]$', fontsize=12)
# plt.tight_layout()
# plt.savefig('ruido_tercios.png')
plt.show()
