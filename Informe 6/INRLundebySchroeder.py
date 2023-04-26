import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
import scipy.signal as sp


# signal = np.sin(np.linspace(0, 2*np.pi, fs))
signal,fs = sf.read('RI-ECM8000-LSS20s.wav', dtype='float32')

"Hilbert transform of a signal"
def hilbert(signal,fs):
    signal=np.abs(np.fft.fft(signal))

    analytic_signal = sp.hilbert(signal)

    amplitude_envelope = np.abs(analytic_signal/np.max(analytic_signal))

    instantaneous_phase = np.unwrap(np.angle(analytic_signal))

    instantaneous_frequency = (np.diff(instantaneous_phase) /(2.0*np.pi) * fs)

    return amplitude_envelope,  instantaneous_frequency, instantaneous_phase

a,b,c=hilbert(signal,fs)
"Plotting the señal con toda la gilada"
# plt.plot(a,np.insert(b,len(b),0))

freqs = np.fft.fftfreq(np.int(len(signal)//2))
print(freqs)
signal=signal/np.max(signal)
signalfft = np.fft.fft(signal[:len(signal)//2])
signalfft=np.real(signalfft)
print(signalfft)
# signal_len=np.int(np.floor(len(signal)/2))
# print(signal_len)
# signalfft=signalfft[signal_len+1:]
# print(len(signalfft))
print(len(freqs))
plt.semilogx(freqs[:((len(freqs)//2)+1)],signalfft[:((len(signalfft)//2)+1)],'r')
# plt.plot(freqs[:((len(freqs)//2)+1)],a[:((len(a)//2))],'b')
plt.show()