import sys
from PyQt5 import uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QMainWindow, QApplication, QFileDialog, QTableWidgetItem, QVBoxLayout, QHBoxLayout
import soundfile as sf
from acoustic_parameters import ac_parameters_mono, ac_parameters_stereo
import sip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas    
import csv

'''
class GUI_Main contains the corresponding programming of the UI

'''

class GUI_Main(QMainWindow):
   
    '''
    the method __init__ starts the objects that will be used for navigate
    through the UI and calculate the acoustical parameters
    
    '''
    
    def __init__(self):
        super().__init__()
        uic.loadUi("UI.ui",self) #This line loads the UI
        self.BotonStereo.setEnabled(False)
        self.BotonMono.clicked.connect(self.leerMONO)
        self.BotonStereo.clicked.connect(self.leerST) #Allows to load the right channel of a STEREO SPLITTED IR
        self.radioRIRMONO.clicked.connect(self.RMono)
        self.radioRIRST.clicked.connect(self.RStereo)
        self.BotonCalc.clicked.connect(self.calcular)
        self.radioRIRBIN.clicked.connect(self.RMono)
        self.comboSuav.currentIndexChanged.connect(self.comboBox)
        if self.comboSuav.currentText()=="Schroeder":
            self.lineVent.setEnabled(False)  
        layout=QHBoxLayout(self.graphicsView) #Defines the plot area location
        self.setLayout(layout)
        self.canvas=FigureCanvas(plt.Figure())
        layout.addWidget(self.canvas) #This lines loads the plot widget canvas
        self.canvas.figure.subplots()
              
    '''
    The method insert_ax loads the graphs axis 
    '''     
    def insert_ax(self):  
        self.ax=self.canvas.figure.subplots()
        self.ax.set_ylim([0, 10])
        self.ax.set_xlim([0, 10])
    '''
    comboBox method enables the window size in case of moving median smoothing 
    '''
    def comboBox(self):
        if self.comboSuav.currentText()=="Schroeder":
            self.lineVent.setEnabled(False)
        if self.comboSuav.currentText()=="Moving Median":
            self.lineVent.setEnabled(True)
    '''
    Method RMono defines the appareance of the buttons that loads the IRs when
    the radio button RIR MONO is clicked
    '''
    def RMono(self):
        self.BotonStereo.setEnabled(False)
        self.BotonMono.setText("Load IR")
    '''
    Method RStereo defines the appareance of the buttons that loads the IRs when
    the radio button RIR STEREO SPLITTED is clicked
    '''    
    def RStereo(self):  
        self.BotonStereo.setEnabled(True)
        self.BotonStereo.setText("Load Right IR")
        self.BotonMono.setText("Load Left IR")

    '''
    Method leerMONO opens a File Dialog that allows to load an IR to calculate
    the Acoustical Parameters in case of clicking the MONO or STEREO IR radio buttons.
    
    '''

    def leerMONO(self):
        archivo=QFileDialog.getOpenFileName(self, 'Open', '.\*') #Opening of a File Dialog
        if self.radioRIRBIN.isChecked() : #In case of loading a STEREO IR, channels Left and Right are loaded
             if archivo[0] !="":
                self.W, self.fs = sf.read(archivo[0])
                self.L=self.W[:,0] #Left Channel
                self.R=self.W[:,1] #Right Channel
                self.W=0
        if self.radioRIRMONO.isChecked() :  #In case of loading a MONO IR         
            self.L=0
            self.R=0
            if archivo[0] !="":
                self.W, self.fs = sf.read(archivo[0])

        if self.radioRIRST.isChecked() :           
            if archivo[0] !="":
                self.L, self.fs = sf.read(archivo[0])
        
        nuevoTexto = archivo[0].split('/')[-1]
        if nuevoTexto:
            self.BotonMono.setText(nuevoTexto)
        
    '''
    Method leerMONO opens a File Dialog that allows to load the right channel of a STEREO SPLITTED IR 
    to calculate the Acoustical Parameters.
    
    '''
    def leerST(self):

        archivo=QFileDialog.getOpenFileName(self, 'Open', '.\*') #Opens a File Dialog to load the right channel of a Stereo Splitted IR
        if self.radioRIRST.isChecked() : 
            self.W=0
            if archivo[0] !="":
                self.R, self.fs = sf.read(archivo[0]) #Reading the right channel of the IR loaded
   
        nuevoTexto = archivo[0].split('/')[-1]
        if nuevoTexto:
            self.BotonStereo.setText(nuevoTexto)
    '''
    Method calcular calls the calculation script "tp10parametros_fix" to obtain the 
    Acoustical Parameters of the corresponding IR loaded
    '''
   
    def calcular(self):
        '''
        Calculation for STEREO IR or STEREO SPLITTED IR
        
        '''
        if np.sum(self.L) != 0 and np.sum(self.R) !=0:
            L=self.L
            R=self.R
            fs=self.fs
            if self.comboFilt.currentText()=="Octave": #In case of applying octave band filter
                ter=0
                nombreColumnas=['31.5 Hz', '63 Hz', '125 Hz', '250 Hz', '500 Hz', '1000 Hz', '2000 Hz', '4000 Hz', '8000 Hz', '16000 Hz'] #Names for the columns of the results table in case of applying octave band filter
                self.tableWidget.setColumnCount(10) #Setting the quantity of columns
                self.tableWidget.setRowCount(8) #Setting the rows quantity
                self.tableWidget.setHorizontalHeaderLabels(nombreColumnas)
            if self.comboFilt.currentText()=="One Third Octave": #In case of applying one third octave band filter
                ter=1
                nombreColumnas=("25 Hz", "31.5 Hz", "40 Hz", "50 Hz", "63 Hz", "80 Hz", "100 Hz", "125 Hz", "160 Hz","200 Hz", "250 Hz", "315 Hz", "400 Hz", " Hz500",
                               "630 Hz", "800 Hz", "1 kHz","1.3 kHz", "1.6 kHz", "2 kHz", "2.5 kHz", "3.2 kHz", "4 kHz", "5 kHz", "6.3 kHz", "8 kHz", "10 kHz", 
                               "12.5 kHz", "16 kHz", "20 kHz") #Names for the columns of the results table in case of applying one third octave band filter
                self.tableWidget.setColumnCount(30) #Setting the quantity of columns
                self.tableWidget.setRowCount(8) #Setting the rows quantity
                self.tableWidget.setHorizontalHeaderLabels(nombreColumnas)
            #The lines 141 to 145 defines the variable "trunc" which corresponds to enable or disable the noise correction
            if self.comboRF.currentText()=="Lundeby":
                trunc=0
            if self.comboRF.currentText()=="Off":
                trunc=1
            #The lines 147 to 150 defines the variable "smooth" which corresponds to enable the chosen smoothing
            if self.comboSuav.currentText()=="Schroeder":
                smooth=0
            if self.comboSuav.currentText()=="Moving Median":
                smooth=1
            nombreFilas=('Tt [s]', 'EDTt [s]', 'C50 [dB]','C80 [dB]','EDT [s]','T20 [s]','T30 [s]','IACCe')
            self.tableWidget.setVerticalHeaderLabels(nombreFilas) #Setting the rows name of the results table
            vent= int(self.lineVent.text()) #Takes the input of the window size chosen in case of Moving Median smoothing
            self.params = ac_parameters_stereo(L, R, ter, trunc, smooth, vent, fs) #Calling the corresponding function for calculate Acoustical Parameters for a STEREO IR
            i=0
            c=np.array(['Tt', 'EDTt', 'C50','C80','EDT','T20','T30','IACCe'])
            columna=0
            for i in range(0,len(c)): #Setting the values into the results table
                for registro in self.params[c[i]]:
                    celda= QTableWidgetItem(str(registro))
                    self.tableWidget.setItem(0, columna, celda)
                    
                    columna += 1
                    
            '''
            Setting the graphs values
            
            '''
                        
            ETC = self.params['ETC'][0]
            smooth = self.params['smooth'][0]
            a = np.arange(0,len(ETC)/fs,1/fs)
            t = a[:len(ETC)]
            ax = self.canvas.figure.axes[0]
            ax.cla()
            ax.plot(t, ETC, label='Energy')
            ax.plot(t, smooth, label=self.comboSuav.currentText())
            try:
                xlim_max = int(np.where(smooth <= -80)[0][0] * 1.1)
                if xlim_max > t.size:
                    xlim_max = t.size-1
            except:
                xlim_max = t.size-1
            ax.set(xlabel='Time [s]', ylabel='Energy [dB]',
                   xlim=(0, t[xlim_max]), ylim=(-100, max(ETC)))
            ax.legend(loc=1)
            ax.figure.tight_layout(pad=0.1)
            self.canvas.draw()
            if self.checkBox.isChecked(): #Exporting the data to a csv file in case of checking the "Export Data" check box
                NombreCsv=QFileDialog.getSaveFileName(self, "Save results","",filter="CSV Files (*.csv)")
                del self.params['ETC']
                del self.params['smooth']
                df=pd.DataFrame(self.params)
                df.index=nombreColumnas
                df.to_csv(NombreCsv[0], sep=';')
            '''
            Calculation for mono IR
            '''
        if np.sum(self.W )!= 0:
            W=self.W
            fs=self.fs
            if self.comboFilt.currentText()=="Octave": #In case of applying octave band filter
                ter=0
                nombreColumnas=['31.5 Hz', '63 Hz', '125 Hz', '250 Hz', '500 Hz', '1000 Hz', '2000 Hz', '4000 Hz', '8000 Hz', '16000 Hz'] #Names for the columns of the results table in case of applying octave band filter
                self.tableWidget.setColumnCount(10) #Setting the quantity of columns
                self.tableWidget.setRowCount(7) #Setting the quantity of rows
                self.tableWidget.setHorizontalHeaderLabels(nombreColumnas)
            if self.comboFilt.currentText()=="One Third Octave": #In case of applying one third octave band filter
                ter=1
                nombreColumnas=("25 Hz", "31.5 Hz", "40 Hz", "50 Hz", "63 Hz", "80 Hz", "100 Hz", "125 Hz", "160 Hz","200 Hz", "250 Hz", "315 Hz", "400 Hz", " Hz500",
                               "630 Hz", "800 Hz", "1 kHz","1.3 kHz", "1.6 kHz", "2 kHz", "2.5 kHz", "3.2 kHz", "4 kHz", "5 kHz", "6.3 kHz", "8 kHz", "10 kHz", 
                               "12.5 kHz", "16 kHz", "20 kHz") #Names for the columns of the results table in case of applying one third octave band filter
                self.tableWidget.setColumnCount(30) #Setting the quantity of columns
                self.tableWidget.setRowCount(7) #Setting the quantity of rows
                self.tableWidget.setHorizontalHeaderLabels(nombreColumnas)
                #The lines 207 to 210 defines the variable "trunc" which corresponds to enable or disable the noise correction
            if self.comboRF.currentText()=="Lundeby":
                trunc=0
            if self.comboRF.currentText()=="Off":
                trunc=1
            #The lines 212 to 215 defines the variable "smooth" which corresponds to enable the chosen smoothing
            if self.comboSuav.currentText()=="Schroeder":
                smooth=0
            if self.comboSuav.currentText()=="Moving Median":
                smooth=1
            nombreFilas=('Tt [s]', 'EDTt [s]', 'C50 [dB]','C80 [dB]','EDT [s]','T20 [s]','T30 [s]')
            self.tableWidget.setVerticalHeaderLabels(nombreFilas) #Setting the rows name of the results table
            vent = int(self.lineVent.text()) #Takes the input of the window size chosen in case of Moving Median smoothing
            self.params = ac_parameters_mono(W, ter, trunc, smooth, vent, fs) #Calling the corresponding function for calculate Acoustical Parameters for a MONO IR
            i = 0
            c = np.array(['Tt', 'EDTt', 'C50','C80','EDT','T20','T30'])
            columna = 0
            for i in range(0, len(c)):#Setting the values into the results table
                for registro in self.params[c[i]]:
                    celda = QTableWidgetItem(str(registro))
                    self.tableWidget.setItem(0, columna, celda)
                    columna += 1
                    
                    
            '''
            Setting the graphs values
            
            '''        
            
            ETC = self.params['ETC']
            smooth = self.params['smooth']
            a = np.arange(0,len(ETC)/fs,1/fs)
            t = a[:len(ETC)]
            ax = self.canvas.figure.axes[0]
            ax.cla()
            ax.plot(t, ETC, label='Energy')
            ax.plot(t, smooth, label=self.comboSuav.currentText())
            try:
                xlim_max = int(np.where(smooth <= -80)[0][0] * 1.1)
                if xlim_max > t.size:
                    xlim_max = t.size-1
            except:
                xlim_max = t.size-1
            ax.set(xlabel='Time [s]', ylabel='Energy [dB]',
                   xlim=(0, t[xlim_max]), ylim=(-100, max(ETC)))
            ax.legend(loc=1)
            ax.figure.tight_layout(pad=0.1)
            self.canvas.draw()
            if self.checkBox.isChecked():
                NombreCsv=QFileDialog.getSaveFileName(self, "Save results","",filter="CSV Files (*.csv)")
                del self.params['ETC']
                del self.params['smooth']
                df=pd.DataFrame(self.params)
                df.index=nombreColumnas
                df.to_csv(NombreCsv[0], sep=';')
'''
Loading the class that calls the UI
'''
if __name__ == '__main__':

    app = QApplication(sys.argv)
    GUI = GUI_Main()
    GUI.show()
    sys.exit(app.exec_())