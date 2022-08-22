# Importación de módulos
from tkinter import *
from Funciones_GUI import *

# Creación de raíz para interfaz gráfica
raiz = Tk()
raiz.title('TP2 - PDS - 1C 2020 - Dragone, Lansky, Lovey')
raiz.minsize(450, 500)
raiz.iconbitmap('untref.ico')

frame1 = LabelFrame(raiz, text='Karplus Strong', padx=5, pady=5)
frame1.pack(padx=5, pady=5)
frame2 = LabelFrame(raiz, text='Finitos delays', padx=5, pady=5)
frame2.pack(padx=5, pady=5)
frame3 = LabelFrame(raiz, text='Infinitos delays', padx=5, pady=5)
frame3.pack(padx=5, pady=5)
frame4 = LabelFrame(raiz, text='Modificaciones de Karplus Strong', padx=5, pady=5)
frame4.pack(padx=5, pady=5)

# Inicio de módulo pygame mixer
mixer.init(44110, 16, 1)

# Creación de botones y cuadros de diálogo, marco 1:
freq_label = Label(frame1, text='Frecuencia en Hz:')
freq_label.grid(row=0, column=0, sticky='e')

frecuencia = IntVar(frame1, value=220)
freq = Entry(frame1, textvariable=frecuencia)
freq.grid(row=0, column=1)

dur_label = Label(frame1, text='Duración en segundos:')
dur_label.grid(row=1, column=0, sticky='e')

duracion = IntVar(frame1, value=2)
dur = Entry(frame1, textvariable=duracion)
dur.grid(row=1, column=1)

fs_label = Label(frame1, text='Frecuencia de muestreo:')
fs_label.grid(row=2, column=0, sticky='e')

sampleRate = IntVar(frame1, value=44100)
samplerate = Entry(frame1, text=sampleRate)
samplerate.grid(row=2, column=1)

KS_bot = Button(frame1, text='Crear y guardar', command=lambda: KS(frecuencia.get(), duracion.get(), sampleRate.get()))
KS_bot.grid(row=0, column=2)

play_bot = Button(frame1, text='Reproducir', command=lambda: play(duracion.get()))
play_bot.grid(row=2, column=2)

# Creación de botones y cuadros de diálogo, marco 2:
nota1 = Label(frame2, text='Ejecutar luego de haber guardado Karplus Strong')
nota1.grid(row=4, column=0, columnspan=2)

fir_guardar = Button(frame2, text='Guardar', command=delay4)
fir_guardar.grid(row=5, column=0)

fir_play = Button(frame2, text='Reproducir', command=lambda: play_fir(duracion.get()))
fir_play.grid(row=5, column=1)

# Creación de botones y cuadros de diálogo, marco 3:
nota1 = Label(frame3, text='Ejecutar luego de haber guardado Karplus Strong')
nota1.grid(row=6, column=0, columnspan=2)

iir_guardar = Button(frame3, text='Guardar', command=delay_inf)
iir_guardar.grid(row=7, column=0)

iir_play = Button(frame3, text='Reproducir', command=lambda: play_iir(duracion.get()))
iir_play.grid(row=7, column=1)

# Creación de botones y cuadros de diálogo, marco 4:
mu = Label(frame4, text='Posición de la cuerda:')
mu.grid(row=8, column=0)

erre = Label(frame4, text='Coeficiente de dinámica:')
erre.grid(row=9, column=0)

barra_mu = Scale(frame4, from_=1, to=99, orient=HORIZONTAL)
barra_mu.grid(row=8, column=1)

barra_r = Scale(frame4, from_=1, to=99, orient=HORIZONTAL)
barra_r.grid(row=9, column=1)

btn_guardar = Button(frame4, text='Guardar', command=lambda: KSM(frecuencia.get(), duracion.get(),
                                                                 barra_mu.get(), barra_r.get(), sampleRate.get()))
btn_guardar.grid(row=8, column=2, rowspan=2)

tuning = Label(frame4, text='Afinación')
tuning.grid(row=10, column=0)

c_tuning = Scale(frame4, from_=1, to=99, orient=HORIZONTAL)
c_tuning.grid(row=10, column=1)

tuning_btn = Button(frame4, text='Guardar afinación', command=lambda: tuning(c_tuning.get()))
tuning_btn.grid(row=10, column=2)

play_btn = Button(frame4, text='Reproducir', command=lambda: play_ksm(duracion.get()))
play_btn.grid(row=11, column=0, columnspan=3, pady=10)


# Ejecución de la ventana
raiz.mainloop()

""" Este programa tiene varios problemas, desde estéticos hasta de control de errores.
Podría llegar a ser una mejor versión, pero el tiempo y la dedicación nos llegaron a un límite.
Gracias por leer hasta el final.

Dragone, Esteban
Lansky, Uriel
Lovey, Milena."""
