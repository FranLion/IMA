#%%
#Imports
from PyQt5.QtWidgets import *
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
from utils import *

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
    def stereocheck(self):
        self.scheck= not self.scheck
        
    def filtercheck(self):
        radioBtn = self.sender()
        if radioBtn.text()=='Tercio':
            self.fcheck = True
        elif radioBtn.text()=='Octavo':
            self.fcheck = False
            
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
        self.sc.draw()
    
    def process_signal(self, signals, stereo=False, iacc=False):
        results = {'t30' : [], 
                    't20' : [],
                    't10' : [],
                    'edt' : [],
                    'C50' : [],
                    'C80' : []}
        for i in range(signals.shape[0]):
            signal = signals[i, :] 
            signal_denoised = denoise(signal, self.fs)
            #Sch
            self.sch_dB = schroeder(signal_denoised)
            params = rt_descriptors(self.sch_dB, signal, self.fs)
            results['t30'].append(params['t30'])
            results['t20'].append(params['t20'])
            results['t10'].append(params['t10'])
            results['edt'].append(params['edt'])
            results['C50'].append(params['C50'])
            results['C80'].append(params['C80'])
            
        if stereo:
            _ = {k+stereo: v for k, v in results.items()}
        
        if isinstance(iacc, np.ndarray):
            results['IACC'] = []
            signals_iacc_L, _ = filtrado(np.squeeze(self.data[:, 0]), self.fs, self.fcheck)
            signals_iacc_R, _ = filtrado(np.squeeze(self.data[:, 1]), self.fs, self.fcheck)
            for i in range(signals_iacc_L.shape[0]):
                signal_L = signals_iacc_L[i, :]
                signal_R = signals_iacc_R[i, :]
                p_iacc = self.IACC_e(signal_L, signal_R, self.fs)
                results['IACC'].append(p_iacc)
        return results
    
    def process_signal_MM(self, signals, stereo=False, iacc=False):
        results = { 't30' : [], 
                    't20' : [],
                    't10' : [],
                    'edt' : [],
                    'C50' : [],
                    'C80' : []}
        for i in range(signals.shape[0]):
            signal = signals[i, :] 
            signal_denoised = denoise(signal, self.fs)
            #Sch
            self.sch_dB=scn.median_filter(self.data,size=50)
            params = rt_descriptors(self.sch_dB, signal, self.fs)
            results['t30'].append(params['t30'])
            results['t20'].append(params['t20'])
            results['t10'].append(params['t10'])
            results['edt'].append(params['edt'])
            results['C50'].append(params['C50'])
            results['C80'].append(params['C80'])
            
        if stereo:
            _ = {k+stereo: v for k, v in results.items()}
        
        if isinstance(iacc, np.ndarray):
            results['IACC'] = []
            signals_iacc_L, _ = filtrado(np.squeeze(self.data[:, 0]), self.fs, self.fcheck)
            signals_iacc_R, _ = filtrado(np.squeeze(self.data[:, 1]), self.fs, self.fcheck)
            for i in range(signals_iacc_L.shape[0]):
                signal_L = signals_iacc_L[i, :]
                signal_R = signals_iacc_R[i, :]
                p_iacc = self.IACC_e(signal_L, signal_R, self.fs)
                results['IACC'].append(p_iacc)
        return results

    def IACC_e(self, L, R, fs):
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
        
        t80 = int(0.08*fs)
        I = np.correlate(L[0:t80], R[0:t80], 'full')/(np.sqrt(np.sum(L[0:t80]*2)*np.sum(R[0:t80]*2)))
        iacce = np.max(np.abs(I)) 
        return iacce
    
    def Calcular(self):       
        
        if self.softcheck=='MMF':
            if self.scheck:
                results = []
                for channel, label in zip(range(self.data.shape[1]), ['L', 'R']):
                    signals, centros = filtrado(np.squeeze(self.data[:, channel]), self.fs, self.fcheck)
                    result = self.process_signal(signals, label, self.data)
                    results.append(result)
                print(results)
            else:
                signals, centros = filtrado(np.squeeze(self.data), self.fs, self.fcheck)
                results = self.process_signal(signals)
        elif self.softcheck=='LS':
            if self.scheck:
                results = []
                for channel, label in zip(range(self.data.shape[1]), ['L', 'R']):
                    signals, centros = filtrado(np.squeeze(self.data[:, channel]), self.fs, self.fcheck)
                    result = self.process_signal(signals, label, self.data)
                    results.append(result)
                print(results)
            else:
                signals, centros = filtrado(np.squeeze(self.data), self.fs, self.fcheck)
                results = self.process_signal(signals)
            
        
        
        print(results)
            # self.sc.fig.clear()
            # axes=self.sc.fig.subplots(1)
            
            # axes.plot(self.sch_dB, label='schroeder')
            # axes.legend()
            # self.sc.draw()
            
        # #Ts
        # self.edt=
        # self.t20=
        # self.t30=



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


def Clarity(filtered_IRs, fs, centros):
    '''
    Calculates the clarity parameters C50 and C80 for a given RIR.
    Parameters
    ----------
    filtered_IRs : array
        Array containing the filtered RIR for each frequency band in each row.
    fs : int
        Sampling frequency.
    Returns
    -------
    C50 : float
    C80 : float
    '''
    
    C50 = np.zeros(centros.size)
    C80 = np.zeros(centros.size)

    
    t50 = np.int64(0.05*fs)  # Index of signal value at 50 ms
    t80 = np.int64(0.08*fs)  # Index of signal value at 80 ms
    
    y50_num = filtered_IRs[0:t50]
    y50_den = filtered_IRs[t50:]
    y80_num = filtered_IRs[0:t80]
    y80_den = filtered_IRs[t80:]
    # y50_num = np.append(y50_num,np.zeros(len(y50_den)-len(y50_num)))
    # y80_num = np.append(y80_num,np.zeros(len(y80_den)-len(y80_num)))

    for banda in range(centros.size):
        c50 = 10*np.log10(np.cumsum(y50_num[banda]) / np.cumsum(y50_den[banda]))
        c80 = 10*np.log10(np.cumsum(y80_num[banda]) / np.cumsum(y80_den[banda]))
        C50[banda]=c50
        C80[banda]=c80

        # c80 = np.round(C50,2)
        # c50 = np.round(C80,2)
        
    return C50, C80

def ttyedtt(x, y, fs):
    '''
    Calculate Tt and early decay time (EDT) from a filtered IR.
    Parameters
    ----------
    x : array
        Array containing the filtered IR per frequency band in each row.
    y : array
        Array containing the .
    fs : int
        Sampling frequency.
    Returns
    -------
    Tt : float
        Transition time.
    EDTt : float
        Early decay time.
    '''
    
    EDTt = []
    Tt = []
    
    for i, ir in enumerate(x):
        
        # Remove the first 5 ms
        
        ir = ir[int(5e-3 * fs):]  
        
        # Find the index of the upper limit of the interval that contains 
        # 99% of the energy 
        
        index = np.where(np.cumsum(ir ** 2) <= 0.99 * np.sum(ir ** 2))[0][-1]
        t_t = index/fs
        Tt.append(t_t)
                
        # Filter the impulse with the moving median filter in order to calculate
        # the parameters
        
        ir2 = y[i]
        ir2 = ir2[:index]  
        t_Tt = np.arange(0, index/fs, 1/fs)
        
        if len(t_Tt) > index:
            t_Tt = t_Tt[:index]
        
        # Calculate minimum squares
            
        A = np.vstack([t_Tt, np.ones(len(t_Tt))]).T
        m, c = np.linalg.lstsq(A, ir2, rcond=-1)[0]
        
        edt_t = -60/m
        EDTt.append(edt_t)
        
    EDTt = np.round(EDTt, 2)
    Tt = np.round(Tt, 2)
    
    return Tt, EDTt

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