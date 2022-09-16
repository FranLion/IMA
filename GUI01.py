#%%
#Imports
import pyqtgraph as pg
from PyQt6.QtWidgets import *
import scipy.signal as sc
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
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes=fig.canvas.draw_idle()
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

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
        texto1=QLabel("Cargar archivos")
        layout.addWidget(texto1)

        boton1=QPushButton("Cargar RI")
        boton1.clicked.connect(self.aguanteFede)
        layout.addWidget(boton1)

        self.botonrad1 = QRadioButton("Octavo")
        self.botonrad1.setChecked(True)
        # self.botonrad1.clicked.connect(filtroBool) 
        # para cuando funcione todo eventualmente hay que 
        # definir un if y seleccionar bandas de octava
        layout.addWidget(self.botonrad1)
            
        self.botonrad2 = QRadioButton("Tercio")
        layout.addWidget(self.botonrad2)
        self.botonrad2.setChecked(True)
            
        boton2=QPushButton("Calcular filtro")
        boton2.clicked.connect(self.aguanteFede2)
        layout.addWidget(boton2)

        sc = MplCanvas(self, width=5, height=4, dpi=100)
        self.sc=sc
        layout.addWidget(self.sc)

# Functions
# %%

    def aguanteFede(self):
        filePath, data, fs, dataMax = loader()
        self.fs=fs
        self.data=data
        self.dataMax=dataMax
        self.sc.axes.cla()  # delete ax1 from the figure
        self.sc.axes.plot(self.data/self.dataMax)
        self.sc.draw()

    def aguanteFede2(self):
        signals = tercio(self.data,self.fs)
        self.signals = signals[0]
        self.signals = np.log10((np.fft.rfft(self.signals)))
        # self.signals = self.signals/np.max(self.data)
        self.mmovil = media_misley(np.array(self.signals),5)
        # self.sc.axes.cla()  # delete ax1 from the figure
        self.sc.axes.plot(self.mmovil)
        self.sc.draw()

## Funciones de cálculos
# %%

def media_misley(x, w):
    """ 
    x: input signal
    w: window
    """
    resultado = np.convolve(x, np.ones(w), 'valid') / w
    return resultado

def loader():
    filePath = QFileDialog.getOpenFileName()
    print(filePath[0])
    data, fs= sf.read(filePath[0])
    data=data.flatten('C') # Ojo aca porque el ruido era stereo asi que hay que cambiar esto segun haga falta
    print(data)
    dataMax=np.max(data)
    indDataMax=np.argmax(data/dataMax)
    data=data[indDataMax:]
    return filePath, data, fs, dataMax

def tercio(data,fs):
    sig = np.fft.fft(data)
    sig = sig[:len(sig)//2]
    indice = np.arange(-16,13,1) # int numbers array [-16,-15,...,12,13]
    fr = 1000 # 
    b = 3
    fm_v = []
    f1_v = []
    f2_v = []
    h_v = []

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
        sos = sc.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',output='sos')
        a,b = sc.butter(6,(f1_v[i],f2_v[i]),btype='bandpass',analog=True,output='ba')
        sig_filtrada = sc.sosfilt(sos, sig)
        signals.append(sig_filtrada)
        # plotwist=plt.plot(signals)
    
    # print(len(signals))
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

# %%
#Closing statements
app = QApplication(sys.argv)
window = Window()
window.show()
sys.exit(app.exec())