# Importación de módulos
from collections import deque
import random as random
import numpy as np
import soundfile as sf
from pygame import mixer
import time


# Karplus Strong original
def KS(freq, dur, fs=44100):
    """
    Sintetiza el sonido de una cuerda de un instrumento como la guitarra y la guarda en .wav.
    :param freq: Frecuencia de la cuerda pulsada.
    :param dur: Duración del audio.
    :param fs: Frecuencia de muestreo, 44100 Hz por defecto.
    :return: Devuevle el vector del audio generado.
    """
    # Datos iniciales
    nSamples = fs * dur
    N = int(fs / freq)

    # Creación de ring buffer
    buf = deque([random.random() - 0.5 for i in range(N)])

    # Creación de vector de ceros
    y = np.array([0] * nSamples, 'float32')

    # Proceso de sintetizador
    for i in range(nSamples):
        y[i] = buf[0]
        avg = 0.5 * (buf[0] + buf[1])
        buf.append(avg)
        buf.popleft()

    # Guardado de audio
    indice = 0
    try:
        for i in range(100):
            with open(f'KS_{i}.wav', 'rb'):
                pass
            indice += 1
            pass
    except FileNotFoundError:
        sf.write(f'KS_{indice}.wav', y, samplerate=fs, subtype='PCM_16')
        pass


# Karplus Strong con LPF y moving pick
def KSM(freq, dur, mu, R, fs=44100):
    """
    Sintetiza el sonido de una cuerda de un instrumento como la guitarra y la guarda en .wav.
    :param freq: Frecuencia de la cuerda pulsada.
    :param dur: Duración del audio.
    :param mu: Distancia del pulsado con respecto al largo de la cueda, entre 0 y 1
    :param R: Coeficiente de LPF y dinámica, entre 0 y 1.
    :param fs: Frecuencia de muestreo, 44100 Hz por defecto.
    :return: Devuevle el vector del audio generado.
    """
    # Datos iniciales
    mu = mu/100
    R = R/100
    nSamples = fs * dur
    N = int(fs / freq)
    aux = int(mu * N)

    # Creación de ring buffer y aplicación de "moving pick"
    buf = np.zeros(N)
    for i in range(N):
        buf[i] = random.random()-0.5
    x1 = np.zeros(N)
    for i in range(N):
        x1[i] = buf[i-aux]
    x1[0:aux] = 0
    aux2 = x1-buf

    # Filtro de dinámica y pasabajos
    aux3 = np.zeros(len(aux2))
    for n in range(len(aux2)):
        aux3[n] = (1-R)*aux2[n]+R*aux2[n-1]
    y1 = deque(aux3)

    # Proceso de sintetizador
    y = np.array([0] * nSamples, 'float32')  # Vector de salida con ceros
    for i in range(nSamples):
        y[i] = y1[0]
        avg = 0.5 * (y1[0] + y1[1])
        y1.append(avg)
        y1.popleft()

    # Guardado de audio
    indice = 0
    try:
        for j in range(100):
            with open(f'KSM_{j}.wav', 'rb'):
                pass
            indice += 1
            pass
    except FileNotFoundError:
        sf.write(f'KSM_{indice}.wav', y, samplerate=fs, subtype='PCM_16')
        pass


# Tuning
def tuning(c):
    """
    Se encarga de afinar los armónicos en caso de frecuencias de muestreo altas
    :param c: Coeficiente de afinación, entre 0 y 1.
    :return: Devuelve un vector del mismo archivo de audio afinado.
    """
    try:
        with open('KSM.wav', 'rb'):
            entrada, fs = sf.read('KSM.wav')
            pass
    except FileNotFoundError:
        print('No se puede abrir el archivo deseado.')
        pass
    c = c/100
    x = np.array(entrada)
    y = np.zeros(len(x))
    for n in range(len(x)):
        y[n] = c*x[n]+x[n-1]-c*y[n-1]

    # Guardado de audio
    indice = 0
    try:
        for i in range(100):
            with open(f'KSM_{i}.wav', 'rb'):
                pass
            indice += 1
            pass
    except FileNotFoundError:
        sf.write(f'KSM_{indice}.wav', y, samplerate=fs, subtype='PCM_16')
        pass


def delay4():
    """Funcion que implementa un numero finito de 4 ecos"""
    x, fs = sf.read('KS_0.wav')
    x = np.array(x)
    alpha = 0.6
    D = 6000
    n = len(x)

    x1 = np.zeros(n)
    x2 = np.zeros(n)
    x3 = np.zeros(n)
    x4 = np.zeros(n)

    for i in range(n):
        x1[i] = x[i-D]
        x2[i] = x[i-2*D]
        x3[i] = x[i-3*D]
        x4[i] = x[i-4*D]

    x1[0:D] = 0
    x2[0:2*D] = 0
    x3[0:3*D] = 0
    x4[0:4*D] = 0

    x1 = x1*alpha
    x2 = x2*alpha**2
    x3 = x3*alpha**3
    x4 = x4*alpha**4

    y = x + x1 + x2 + x3 + x4
    y = y / max(abs(y))
    # Guardado de audio
    indice = 0
    try:
        for i in range(100):
            with open(f'KS_FIR_{i}.wav', 'rb'):
                pass
            indice += 1
            pass
    except FileNotFoundError:
        sf.write(f'KS_FIR_{indice}.wav', y, samplerate=fs, subtype='PCM_16')
        pass


def delay_inf():
    """Funcion que implementa un numero infinito de ecos"""
    x, fs = sf.read('KS_0.wav')
    x = np.array(x)
    alpha = 0.7
    D = 6000
    n = len(x)
    y = np.zeros(n)

    for i in range(n):
        y[i] = x[i] + alpha * y[i - D]

    y = y / max(abs(y))
    # Guardado de audio
    indice = 0
    try:
        for i in range(100):
            with open(f'KS_IIR_{i}.wav', 'rb'):
                pass
            indice += 1
            pass
    except FileNotFoundError:
        sf.write(f'KS_IIR_{indice}.wav', y, samplerate=fs, subtype='PCM_16')
        pass


def play(tiempo):
    mixer.init()
    mixer.music.load('KS_0.wav')
    mixer.music.play()
    time.sleep(tiempo)
    mixer.quit()


def play_fir(tiempo):
    mixer.init()
    mixer.music.load('KS_FIR_0.wav')
    mixer.music.play()
    time.sleep(tiempo)
    mixer.quit()


def play_iir(tiempo):
    mixer.init()
    mixer.music.load('KS_IIR_0.wav')
    mixer.music.play()
    time.sleep(tiempo)
    mixer.quit()


def play_ksm(tiempo):
    mixer.init()
    mixer.music.load('KSM_0.wav')
    mixer.music.play()
    time.sleep(tiempo)
    mixer.quit()
