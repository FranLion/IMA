#%%
#Imports
from PyQt6.QtWidgets import *
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

## clases GUI
# %%
#Canvas class
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes=self.fig.canvas.draw_idle()
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

#Window class
class Window(QWidget):
    def __init__(self):
        super(Window,self).__init__()
        # self.setWindowIcon(QIcon("Icono.jpg"))
        self.setWindowTitle("IMA - Instrumentos y mediciones Acústicas - Adquisición de RIR ")
        self.resize(1200,700)

        #Layout
        layout=QVBoxLayout()
        self.setLayout(layout)
        #Created objects
        layout_uptexts=QGridLayout()
        texto1=QLabel("Cargar archivos")
        layout_uptexts.addWidget(texto1)
        self.texto2=QLabel("[El archivo cargado]")
        layout_uptexts.addWidget(self.texto2,0,3)
        layout.addLayout(layout_uptexts)


        loader_layout = QHBoxLayout()
        layout.addLayout(loader_layout)
        self.botonRIMono=QPushButton("Cargar RI Mono")
        self.botonRIMono.clicked.connect(self.cargarRIs)
        loader_layout.addWidget(self.botonRIMono)
        
        pholder = QHBoxLayout()
        layout.addLayout(pholder)

        self.scheck=False
        self.botonRIStereo=QPushButton("Cargar RI Stereo")
        self.botonRIStereo.clicked.connect(self.stereocheck)
        self.botonRIStereo.clicked.connect(self.cargarRIs)

        loader_layout.addWidget(self.botonRIStereo)

        layout1=QHBoxLayout()
        layout.addLayout(layout1)

        filter_groupbox = QGroupBox('Filtering')
        filter_box = QVBoxLayout()


        self.botonrad1 = QRadioButton("Octavo")
        self.botonrad1.toggled.connect(self.filtercheck)
        self.botonrad1.setChecked(True)
        self.botonrad2 = QRadioButton("Tercio")
        self.botonrad2.toggled.connect(self.filtercheck)
        # self.botonrad2.setChecked(True)

        filter_groupbox2 = QGroupBox('Smoothing')
        filter_box2 = QVBoxLayout()

        self.botonrad3 = QRadioButton("Lundeby + Schroeder")
        self.botonrad3.toggled.connect(self.smoothcheck)
        self.botonrad3.setChecked(True)
        self.botonrad4 = QRadioButton("MMF")
        self.botonrad4.toggled.connect(self.smoothcheck)
        # self.botonrad4.setChecked(True)

        filter_box2.addWidget(self.botonrad3)
        filter_box2.addWidget(self.botonrad4)
        filter_groupbox2.setLayout(filter_box2)
        
        layout1.addWidget(filter_groupbox)
        layout1.addWidget(filter_groupbox2)
        
        filter_box.addWidget(self.botonrad1)
        filter_box.addWidget(self.botonrad2)
        filter_groupbox.setLayout(filter_box)

        # boton3 = QRadioButton('un culo que ver')
        # layout.addWidget(boton3)

        self.data = None    
        boton2=QPushButton("Calcular")
        boton2.clicked.connect(self.Calcular)
        layout.addWidget(boton2)

        sc = MplCanvas(self, width=5, height=3, dpi=100)
        self.sc=sc
        layout.addWidget(self.sc)
        
        tablaca=QTableWidget()
        tablaca.setHorizontalHeaderLabels(['Tt [s]', 'Tt [s]', 'Tt [s]', 'Tt [s]', 'Tt [s]', 'Tt [s]'])
        tablaca.show
        layout.addWidget(tablaca)

# Functions
# %%
    def stereocheck(self):
        self.scheck= not self.scheck
    def filtercheck(self):
        radioBtn = self.sender()
        if radioBtn.text()=='Tercio':
            bands = [ 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,10000, 12500, 16000]
            # bands = [10, 12 , 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000, 25000, 31500,40000,50000,63000,80000,100000,125000,160000,200000]
            self.bands=np.array(bands)
            print(self.bands)
        elif radioBtn.text()=='Octavo':
            bands = [16,31.5,63,125,250,500,1000,2000,4000,8000,16000]
            self.bands= np.array(bands)
            print(self.bands)    
    def smoothcheck(self):
        radioBtn = self.sender()
        if radioBtn.text()=='MMF':
            self.softcheck='MMF'
        elif radioBtn.text()=='Lundeby + Schroeder':
            self.softcheck='LS'

    def cargarRIs(self):
        file_name, data, fs, data4filter = loader(self.scheck)
        print(self.scheck)
        print(file_name)
        if not file_name:
            return
        self.texto2.setText(file_name + "  sucessfully loaded")
        self.fs=fs
        self.data=data
        self.file_name=file_name
        self.sc.fig.clear()
        
        if self.scheck==True:
            number_channels = data.shape[1]
            axes = self.sc.fig.subplots(nrows=number_channels, squeeze=False)
            for axis, row in zip(axes, data.T):
                axis[0].plot(row)
        else:
            number_channels = data.shape[0]
            axes=self.sc.fig.subplots(1)
            axes.plot(self.data)
            # plotdata=np.array(self.data)
            # self.sc.axes.plot(self.data)

        self.scheck=False
        self.sc.draw()

    def Calcular(self):
        
        signals = filtrado(self.file_name,self.bands)

        # print(self.mmovil)
        # for axis,row in zip(self.sc.fig.axes,self.signals.T):
        #     axis.plot(row)

        if self.softcheck=='MMF':
            test=scn.median_filter(signals,size=9)
            print("señal filtrada: ",test)
        elif self.softcheck=='LS':

            # etc = E_norm(signals)
            t_env,etc=envelope(signals, self.fs)
            print("etc",etc)
            # signals = 10 * np.log10(ETC + sys.float_info.epsilon)
            etc=np.array(etc)
            #Lu
            # cruce,c=lundeby(etc,self.fs)
            # cruce = 10*np.log10(cruce)
            # print("Cruce: ", cruce)
            #Sch
            self.sch_dB = schroeder(signals)
            # print(self.sch_dB)
            self.sc.fig.clear()
            axes=self.sc.fig.subplots(1)
            # axes.plot(cpoint)
            # axes.plot(t_env,self.env, label='env')
            
            axes.plot(self.sch_dB)
            # axes.plot(self.sch_dB[0:cruce], label='schro')
            # x=np.arange(0,len(self.sch_dB))
            # axes.axhline(cruce, label='Loco Lundeby',linestyle='-', color = 'r',)
            axes.legend()
            # axes.plot(t_env,cruce)
            self.sc.draw()
            # self.sch_dB=np.nan_to_num(self.sch_dB, copy=True, nan=0.0, posinf=None, neginf=None)
            # print("El pichichi Schroeder: ", self.sch_dB)
            # print("limite superior: ",limsup)

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
    x.flatten()
    ETC = np.zeros(x.shape)  
    
    for i, y in enumerate(x):
        E = np.abs(sc.hilbert(y))**2
        ETC[i] = E/np.max(E)
        
    return ETC

def schroeder(IR):
    # Schroeder integration
    sch = np.cumsum(IR[::-1]**2)[::-1]
    sch_dB = 10.0 * np.log10(sch / np.max(sch))
    return sch_dB

def media_misley(x, w):
    """ 
    x: input signal
    w: window
    """
    if x.shape[1] > 1:
        result = None
        for array in x.T:
            column = np.column_stack([array])
            filtered=media_misley(column, w)
            if result is None:
                result = np.empty((filtered.shape[0], 0))
            result=np.append(result,filtered,axis=1)
        return result
    resultado = np.convolve(x.flatten(), np.ones(w), 'full') / w
    return resultado[:,np.newaxis]

def loader(boolCH):
    if not boolCH:
        boolCH=False
    file_name, _ = QFileDialog.getOpenFileName()
    if not file_name:
        return [None] * 4
    data, fs= sf.read(file_name,always_2d=True)
    data4filter=[data,fs]
    # dataMax=np.max(abs(data))
    # data=data/dataMax
    indDataMax=np.argmax(data)
    data=data[indDataMax:]
    return file_name, data, fs, data4filter

def filtrado(file_name, bands):  # pylint: disable=too-many-locals
    raw_signal,fs = sf.read(file_name,always_2d=True)
    band_type = ac._check_band_type(bands)

    if band_type == 'octave':
        low = ac.octave_low(bands[0], bands[-1])
        high = ac.octave_high(bands[0], bands[-1])
        maxband=np.max(low)
    elif band_type == 'third':
        low = ac.third_low(bands[0], bands[-1])
        high = ac.third_high(bands[0], bands[-1])
        maxband=1
    
    print(maxband)
    for band in range(bands.size):
        # Filtering signal
        filtered_signal = ac.bandpass(raw_signal, low[band]/maxband, high[band]/maxband, fs, order=8)
        abs_signal = np.abs(filtered_signal) / np.max(np.abs(filtered_signal))
        # fitlered_signal=
    # print(abs_signal)
    return abs_signal


def estim_slope(t_env, env, init, end):
    # Solo me interesa el valor de la slope
    init_idx = np.where(env < init)[0][0]
    try:
        end_idx = np.where(env < end)[0][0]
    except:
        end_idx = len(env)-1
    # regresion lineal 
    x = t_env[init_idx:end_idx+1]
    y = env[init_idx:end_idx+1]
    slope, intercept = st.linregress(x,y)[0:2]
    return intercept, slope , x, y

# def lundeby(signal, fs, time_interval=0.0050):
#     # CONSTANTES DE DISEÑO
#     TIME_INTERVAL = 0.009 # [s] from 0.005 to 0.001 
#     NOISE_FLOOR_DISTANCE = 5 # [dB] from 5 to 10. Level above the noise
#     INTERVALS = 10 # Intervals per 10 dB of decay. From 3 to 10 for low - high freqs
#     MARGIN = 5 # Safety margin from cross point. From 5 to 10 dB of decay 
#     DINAMIC_ABOVE, DINAMIC_BELOW = 10, 5 # Dinamic range of 10-20 dB referred to the noise floor
    
#     interval = int(time_interval * fs) 
    
#     n_windows = len(signal) // interval 
#     remainder = len(signal) % interval
#     env = np.empty(n_windows)
#     for i in range(n_windows):
#         env[i] = signal[i*interval:(i+1)*interval].sum()/interval
#     env = env / np.max(abs(env))
#     t_env = np.arange(0,n_windows*interval, interval)
#     env = amplitude_to_db(env)
#     env = np.array(env, dtype='int32')
#     # standarization
#     # onset = np.argmax(abs(signal))
#     # signal = signal[onset:]
#     # signal = signal / np.max(abs(signal))
#     init=signal[1]
#     end=signal[-1]
    
#     # squared response
#     #signal_sqr = np.power(signal, 2)
#     signal_sqr = abs(signal)
#     t = np.arange(0,len(signal_sqr))

#     # average smoothing
#     # t_env, env = envelope(signal_sqr, fs, time_interval=TIME_INTERVAL)
    
#     # First estimation of noise floor using the tail (last 10%)
#     tail = int(len(t_env) * 0.1)
#     noise_level = env[-tail:].sum() / tail
#     #print('First estimation of noise floor: {:.2f} dB'.format(noise_level))

#     # intercept, slope, x_line, y_line = estim_slope(t_env, env, 0, noise_level+NOISE_FLOOR_DISTANCE)
    
#     init_idx = np.where(env < init)
#     try:
#         end_idx = np.where(env < end)
#     except:
#         end_idx = len(env)-1
#     # regresion lineal 
#     x = t_env[init_idx:end_idx+1]
#     y = env[init_idx:end_idx+1]
#     slope, intercept = st.linregress(x,y)[0:2]

#     cross_point = (noise_level - intercept) / slope

#     # Find new time interval
#     intervals_per_10dB = 6 #3 - 10 [low - high]
#     interval_dB = 10 / intervals_per_10dB
#     interval = np.int32(-interval_dB / slope)
#     time_interval = interval / fs
#     #print('New time interval: {:.4f} seconds'.format(time_interval))

#     t_env, env= envelope(signal_sqr, fs, time_interval=time_interval)
    
#     for i in range(5):
#         margin_cross = 7 #5-10dB
#         safe_cross_point = int(-margin_cross/slope) + int(cross_point)
#         tail = int(len(t_env) * 0.1)
#         if (safe_cross_point < t_env[-tail]):
#             #print('uso el intervalo')
#             index_cross = np.where(t_env > safe_cross_point)[0]
#             noise_level = env[index_cross:].sum() / len(env[index_cross:])
#         else:
#             #print('uso la tail')
#             noise_level = env[-tail:].sum() / tail
#         #print('Nueva estimacion del piso de ruido de {:.2f} dB'.format(noise_level))


#         def estim_slope_f(t_env, env, init, end):
#             x = t_env[init:end+1]
#             y = env[init:end+1]
#             slope, intercept = sc.stats.linregress(x,y)[0:2]
#             return intercept, slope , x, y


#         # Estimar la pendiente 5 dB [5-10] encima del piso de ruido para un rango de 10 dB [10-20]

#         init = (noise_level + 10 - intercept) / slope
#         if init < 0 :
#             init = 0 
#         init = int(init / (time_interval * fs))
#         end = (noise_level-5 - intercept) / slope
#         end = int(end / (time_interval * fs))

#         intercept_f, slope_f, x_line_f, y_line_f = estim_slope_f(t_env, env, init, end)
#         cross_point = (noise_level - intercept_f) / slope_f
#     # insert delay samples
#     cross_point = cross_point + init
#     return int(cross_point)

def envelope(signal, fs, time_interval=0.0050):
    #hacer time interval variable por banda
    #time_interval = 0.0050 #10 - 50 ms
    interval = int(time_interval * fs) 
    
    n_windows = len(signal) // interval 
    remainder = len(signal) % interval
    env = np.empty(n_windows)
    for i in range(n_windows):
        env[i] = signal[i*interval:(i+1)*interval].sum()/interval
    env = env / np.max(abs(env))
    t_env = np.arange(0,n_windows*interval, interval)
    return t_env, amplitude_to_db(env)

def lundeby(IR, Fs):
   
    N = IR.size
    energy = IR
    # energy=energy()
    med = np.zeros(np.int32(N/(Fs*0.01)),dtype='int32')
    eje_tiempo = np.zeros(np.int32(N/(Fs*0.01)),dtype='int32')
    
# Divide in sections and calculate the mean.    
    t = np.floor(N/(Fs*0.01)).astype('int')
    v = np.floor(N/t).astype('int')   
    for i in range(0, t):
        med[i] = np.mean(energy[i * v:(i + 1) * v])
        eje_tiempo[i] = np.ceil(v/2).astype('int') + (i*v)
        

# Calculate noise level of the last 10% of the signal.    
    rms_dB = 10 * np.log10(np.sum(energy[np.int32(np.round(0.9 * N)):]) / (0.1 * N) / np.max(energy))
    meddB = 10 * np.log10(med / np.max(energy))

# The linear regression of the 0dB interval and the mean closest to the noise + 10dB is sought.   
    try:
        r = int(np.max(np.argwhere(meddB > rms_dB + 10)))
           
        if np.any(meddB[0:r] < rms_dB+10):
            r = np.min(np.min(np.where(meddB[0:r] < rms_dB + 10)))
        if np.all(r==0) or r<10:
            r=10
    except:
        r = 10

# Least squares.       
    A = np.vstack([eje_tiempo[0:r], np.ones(len(eje_tiempo[0:r]))]).T
    m, c = np.linalg.lstsq(A, meddB[0:r], rcond=-1)[0]
    cruce = np.int32((rms_dB-c)/m)
    
# Insufficient SNR.    
    if rms_dB > -20:        
        punto = len(energy)
        C = None        
    else:
        error = 1
        INTMAX = 50
        veces = 1               
        while error > 0.0004 and veces <= INTMAX:
            
# Calculates new time intervals for the mean with approximately p steps for each 10 dB.            
            p = 10
            
# Number of samples for the decay slope of 10 dB.            
            delta = np.int(abs(10/m)) 
            
# Interval over which the mean is calculated.           
            v = np.floor(delta/p).astype('int') 
            t = int(np.floor(len(energy[:int(cruce-delta)])/v))            
            if t < 2:
                t = 2
                
            elif np.all(t == 0):
                t = 2
            media = np.zeros(t)
            eje_tiempo = np.zeros(t)
            
            for i in range(0, t):
                media[i] = np.mean(energy[i*v:(i + 1) * v])
                eje_tiempo[i] = np.ceil(v / 2) + (i * v).astype('int')
                
            mediadB = 10 * np.log10(media / max(energy))
            A = np.vstack([eje_tiempo, np.ones(len(eje_tiempo))]).T
            m, c = np.linalg.lstsq(A, mediadB, rcond=-1)[0]

# New noise average level, starting from the point of the decay curve, 10 dB below the intersection.           
            noise = energy[int(abs(cruce + delta)):]
            
            if len(noise) < round(0.1 * len(energy)):
                noise = energy[round(0.9 * len(energy)):]
                
            rms_dB = 10 * np.log10(sum(noise)/ len(noise) / np.max(energy))

# New intersection index           
            error = abs(cruce - (rms_dB - c) / m) / cruce
            cruce = np.round((rms_dB - c) / m)
            veces += 1
                   
# Output validation
    if cruce > N:
        punto = N
    else:
        punto = int(cruce)       
    C = np.max(energy) * 10 ** (c / 10) * np.exp(m/10/np.log10(np.exp(1))*cruce) / (
        -m / 10 / np.log10(np.exp(1)))
        
    return punto, C

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
        slope, intercept = ac.room.stats.linregress(x, y)[0:2]

        # Reverberation time (T30, T20, T10 or EDT)
        db_regress_init = (init - intercept) / slope
        db_regress_end = (end - intercept) / slope
        t60[band] = factor * (db_regress_end - db_regress_init)
    return t60

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

#Closing statements
app = QApplication(sys.argv)
window = Window()
window.show()
sys.exit(app.exec())