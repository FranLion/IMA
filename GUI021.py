#%%
#Imports
from ctypes import alignment
import pyqtgraph as pg
from PyQt6.QtWidgets import *
import scipy.signal as sc
import scipy.ndimage as scn
import acoustics
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
        self.resize(1200,600)

        #Layout
        layout=QVBoxLayout()
        self.setLayout(layout)
        #Created objects
        layout_uptexts=QGridLayout()
        texto1=QLabel("Cargar archivos")
        layout_uptexts.addWidget(texto1)
        texto2=QLabel("[El archivo cargado]")
        layout_uptexts.addWidget(texto2,0,3)
        layout.addLayout(layout_uptexts)


        loader_layout = QHBoxLayout()
        layout.addLayout(loader_layout)
        self.botonRIMono=QPushButton("Cargar RI Mono")
        self.botonRIMono.clicked.connect(self.cargarRIs)
        loader_layout.addWidget(self.botonRIMono)


        self.scheck=False
        self.botonRIStereo=QPushButton("Cargar RI Stereo")
        self.botonRIStereo.clicked.connect(self.stereocheck)
        self.botonRIStereo.clicked.connect(self.cargarRIs)

        loader_layout.addWidget(self.botonRIStereo)

        layout1=QHBoxLayout()
        layout.addLayout(layout1)

        filter_groupbox = QGroupBox('Filtering')
        filter_box = QVBoxLayout()

        botonrad1 = QRadioButton("Octavo")
        botonrad2 = QRadioButton("Tercio")
        botonrad2.setChecked(True)

        filter_groupbox2 = QGroupBox('Smoothing')
        filter_box2 = QVBoxLayout()

        botonrad3 = QRadioButton("Lundeby + Schroeder")
        botonrad4 = QRadioButton("MMF")
        botonrad4.setChecked(True)

        filter_box2.addWidget(botonrad3)
        filter_box2.addWidget(botonrad4)
        filter_groupbox2.setLayout(filter_box2)
        layout1.addWidget(filter_groupbox2)
        
        filter_box.addWidget(botonrad1)
        filter_box.addWidget(botonrad2)
        filter_groupbox.setLayout(filter_box)
        layout1.addWidget(filter_groupbox)

        boton3 = QRadioButton('un culo que ver')
        layout.addWidget(boton3)

        self.data = None    
        boton2=QPushButton("Calcular filtro")
        boton2.clicked.connect(self.calculoFiltros)
        layout.addWidget(boton2)

        sc = MplCanvas(self, width=5, height=4, dpi=100)
        self.sc=sc
        layout.addWidget(self.sc)
        
        tablaca=QTableWidget()
        layout.addWidget(tablaca)

# Functions
# %%
    def stereocheck(self):
        self.scheck= not self.scheck

    def cargarRIs(self):
        filePath, data, fs, dataMax = loader(self.scheck)
        print(self.scheck)
        if not filePath:
            return
        self.fs=fs
        self.data=data
        self.dataMax=dataMax
        self.sc.fig.clear()
        
        # frec=np.fft.fft(self.data)
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

    def calculoFiltros(self):
        if self.data is None:
            return
        self.signals = tercio(self.data,self.fs)
        self.signals = np.array(self.signals)
        # for i, signal in enumerate(signals):
        #     sf.write(f"señal{i}.wav",signal,self.fs)
        self.mmovil = scn.median_filter(self.signals,size=3)

        print(self.mmovil)
        for axis,row in zip(self.sc.fig.axes,self.signals.T):
            axis.plot(row)

        self.sc.draw()

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
    filePath, _ = QFileDialog.getOpenFileName()
    if not filePath:
        return [None] * 4
    data, fs= sf.read(filePath,always_2d=True)
    dataMax=np.max(abs(data))
    # data=data/dataMax
    indDataMax=np.argmax(data)
    data=data[indDataMax:]
    return filePath, data, fs, dataMax

def tercio(data,fs):
    indice = np.arange(-16,13,1) # int numbers array [-16,-15,...,12,13]
    fr = 1000 # 
    b = 3
    fm_v = []
    f1_v = []
    f2_v = []

    signals = []

    # if data.shape[1]>1:
    #     result = None
    #     for array in data

    for i in range(0,len(indice)):
        fm = fr*2**(indice[i]/b)
        f1 = fm*(2**(-1/(2*b)))/(0.5*fs)
        f2 = fm*(2**(1/(2*b)))/(0.5*fs)
        fm_v.append(fm)
        f1_v.append(f1)
        f2_v.append(f2)

    for i in range(0,len(indice)):
        sos = sc.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',output='sos')
        a,b = sc.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',analog=True,output='ba')
        sig_filtrada = sc.sosfilt(sos, data)
        signals.append(sig_filtrada)
    
    return signals

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