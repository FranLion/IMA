import numpy as np
import soundfile as sf

RIR, fs = sf.read("")
RIR_max = np.where(RIR == np.max(RIR))[0][0]
RIR = RIR[int(RIR_max - (0.5 * fs)): int(RIR_max + (3 * fs)):]
RMS = lambda signal: np.sqrt(np.mean(signal**2))
dB = lambda signal: 20 * np.log10(abs(signal))

# Se define el nivel de ruido a partir del último 10% de la señal
noise = RIR[round(0.9 * len(RIR)):len(RIR)]
noise_RMS = RMS(noise)
LN = dB(noise_RMS)
LIR = dB(np.max(abs(RIR)))
INR = LIR - LN