import matplotlib.pyplot
import numpy as np
import scipy.io.wavfile
import scipy.signal
import scipy.fftpack

fs = 48000

T = 4

t = np.arange(0, int(T*fs)) / fs

log_signal = scipy.signal.chirp(t, f0=1, f1=20000, t1=T, method='logarithmic')
loginv_signal=log_signal[::-1]


log_n = len(log_signal) # length of the signal
log_k = np.arange(log_n)
log_T = log_n/fs
log_frq = log_k/log_T # two sides frequency range
log_frq = log_frq[range(log_n//2)] # one side frequency range
log_Y = np.fft.fft(log_signal)/log_n # fft computing and normalization
log_Y = log_Y[range(log_n//2)]

loginv_n = len(loginv_signal) # length of the signal
loginv_k = np.arange(loginv_n)
loginv_T = loginv_n/fs
loginv_frq = loginv_k/loginv_T # two sides frequency range
loginv_frq = loginv_frq[range(loginv_n//2)] # one side frequency range
loginv_Y = np.fft.fft(loginv_signal)/loginv_n # fft computing and normalization
loginv_Y = loginv_Y[range(loginv_n//2)]


fig, axes = matplotlib.pyplot.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,figsize=(10,5))
axes[0].set_title('Logarithmic Sweep')
axes[0].plot(log_signal)
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Amplitude')
axes[0].set_ylim(-1,1)
axes[1].set_title('Inverse Logarithmic Sweep')
axes[1].plot(loginv_signal)
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Amplitude')
axes[1].set_ylim(-1,1)
matplotlib.pyplot.savefig('time_domain_compare_logs_New.png', bbox_inches="tight")
matplotlib.pyplot.show()

# convolve
impulse_response = scipy.signal.fftconvolve(log_signal, loginv_signal, mode='same')

IR_n = len(impulse_response) # length of the signal
IR_k = np.arange(IR_n)
IR_T = IR_n/fs
IR_frq = IR_k/IR_T # two sides frequency range
IR_frq = IR_frq[range(IR_n//2)] # one side frequency range
IR_Y = np.fft.fft(impulse_response)/IR_n # fft computing and normalization
IR_Y = IR_Y[range(IR_n//2)]

fig, axes = matplotlib.pyplot.subplots(2, 1, sharex=False, sharey=False, constrained_layout=True,figsize=(10,5))
axes[0].plot(impulse_response,'r') # plotting the spectrum
axes[0].set_title('Time Domain')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Amplitude')
axes[0].set_xlim(len(impulse_response)*(0.5-0.005),len(impulse_response)*(0.5+0.005))
axes[1].semilogx(log_frq, 20*np.log10(abs(IR_Y)),'r') # plotting the spectrum
axes[1].set_title('Frequency Domain')
axes[1].set_xlabel('Frequency [Hz]')
axes[1].set_ylabel('Magnitude [dB]')
axes[1].set_xlim(20,20000)
axes[1].set_ylim(-120,0)
matplotlib.pyplot.savefig('response_of_IR_New.png', bbox_inches="tight")
matplotlib.pyplot.show()