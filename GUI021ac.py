#%%
#Imports
from ctypes import alignment
import pyqtgraph as pg
from PyQt6.QtWidgets import *
import scipy.signal as sc
import scipy.ndimage as scn
import acoustics.room as ac
import soundfile as sf
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
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
        self.botonrad2 = QRadioButton("Tercio")
        self.botonrad2.toggled.connect(self.filtercheck)
        # self.botonrad2.setChecked(True)

        filter_groupbox2 = QGroupBox('Smoothing')
        filter_box2 = QVBoxLayout()

        self.botonrad3 = QRadioButton("Lundeby + Schroeder")
        self.botonrad3.toggled.connect(self.smoothcheck)
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
        # print(signals)
        
        # self.signals = np.array(self.signals)
        # # for i, signal in enumerate(signals):
        # #     sf.write(f"señal{i}.wav",signal,self.fs)
        # self.mmovil = scn.median_filter(self.signals,size=3)

        # print(self.mmovil)
        # for axis,row in zip(self.sc.fig.axes,self.signals.T):
        #     axis.plot(row)

        if self.softcheck=='MMF':
            test=scn.median_filter(signals,size=9)
            print("señal filtrada: ",test)
        elif self.softcheck=='LS':
            #Sch
            self.sch_dB = schroeder(signals)
            print(self.sch_dB)
            axes=self.sc.fig.subplots(1)
            axes.plot(self.sch_dB)
            #Lu
            punto,C=lundeby(self.sch_dB,self.fs)
            print("Limite superior: ",punto)
            print("Cruce: ",C)

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
    
    # print(abs_signal)
    return abs_signal

def lundeby(IR, Fs):
    IR=IR[1:]
    N = IR.size
    energy = IR
    med = np.zeros(int(N/(Fs*0.01)))
    eje_tiempo = np.zeros(int(N/(Fs*0.01)))
    enmax=np.max(energy)

# Divide in sections and calculate the mean.    
    t = np.floor(N/(Fs*0.01)).astype('int')
    v = np.floor(N/t).astype('int')   
    for i in range(0, t):
        med[i] = np.mean(energy[i * v:(i + 1) * v])
        eje_tiempo[i] = np.ceil(v/2).astype('int') + (i*v)
        
# Calculate noise level of the last 10% of the signal.    
    rms_dB = 10 * np.log10(np.sum(energy[round(0.9 * N):]) / (0.1 * N) / enmax)
    meddB = 10 * np.log10(med / enmax)

# The linear regression of the 0dB interval and the mean closest to the noise + 10dB is sought.   
    try:
        r = int(max(np.argwhere(meddB > rms_dB + 10)))
           
        if np.any(meddB[0:r] < rms_dB+10):
            r = min(min(np.where(meddB[0:r] < rms_dB + 10)))
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
            delta = int(abs(10/m)) 
            
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
                
            mediadB = 10 * np.log10(media / enmax)
            A = np.vstack([eje_tiempo, np.ones(len(eje_tiempo))]).T
            m, c = np.linalg.lstsq(A, mediadB, rcond=-1)[0]

# New noise average level, starting from the point of the decay curve, 10 dB below the intersection.           
            noise = energy[int(abs(cruce + delta)):]
            
            if len(noise) < round(0.1 * len(energy)):
                noise = energy[round(0.9 * len(energy)):]
                
            rms_dB = 10 * np.log10(sum(noise)/ len(noise) / enmax)

# New intersection index           
            error = abs(cruce - (rms_dB - c) / m) / cruce
            cruce = round((rms_dB - c) / m)
            veces += 1

# Output validation            
    if cruce > N:
        punto = N
    else:
        punto = int(cruce)       
    C = enmax * 10 ** (c / 10) * np.exp(m/10/np.log10(np.exp(1))*cruce) / (
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