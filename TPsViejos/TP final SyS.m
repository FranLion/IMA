%% VENTANA INICIO

function New_Inicio(~,~)
close all

%% Crear una figura para usar la GUI

f1 = figure('Visible','off','Position',[0 0 300 410],'Resize','off');
    movegui(f1,'center');
    f1.Name = 'Inicio';
    f1.MenuBar = 'none';
    f1.NumberTitle = 'off';
    f1.Visible = 'on';

%% Crear paneles

P1 = uipanel (f1,'title', '',...
                 'units','normalized','Position',[.068 .046 .858 .92]); 
    
%% Crear pregunta

handles.NI.txt1 = uicontrol(P1,'style','text','FontSize',14,...
                  'string','¿Qué desea hacer?',...
                  'units','normalized','position',[.08 .67 .848 .183]);
              
%% Crear botones de interacción

handles.btnRR = uicontrol('parent',P1,...
                  'style','pushbutton',...
                  'string','Ajustar el dodecaedro',...
                  'units','normalized','position',[.08 .5 .848 .094],...
                  'callback',@RRajuste);

handles.btnSS = uicontrol('parent',P1,...
                  'style','togglebutton',...
                  'string','Medir el TR',...
                  'units','normalized','position',[.08 .35 .848 .094],...
                  'callback',@SSmedicion);

handles.btnTR = uicontrol('parent',P1,...
                  'style','pushbutton',...
                  'string','Calcular el TR',...
                  'units','normalized','position',[.08 .2 .848 .094],...
                  'callback',@CalculoTR);   
             
              
end

%% VENTANA RUIDO ROSA

function RRajuste(~,~)
close all; 

%% Crear una figura para usar la GUI

f3 = figure('Visible','off','Position',[0 0 300 410],'Resize','off');
    movegui(f3,'center');
    f3.Name = 'Ajuste';
    f3.MenuBar = 'none';
    f3.NumberTitle = 'off';
    f3.Visible = 'on'; 
        
%% Crear paneles

rrP1 = uipanel (f3,'title', '',...
                 'units','normalized','Position',[.068 .363 .858 .586]);
rrP2 = uipanel (f3,'title', '',...
                 'units','normalized','Position',[.068 .046 .858 .297]);
             
%% Crear títulos

handles.RR.txt1 = uicontrol(rrP1,'style','text','FontSize',14,...
                  'string','Ajuste de la fuente con ruido rosa',...
                  'units','normalized','position',[.08 .71 .848 .183]);

handles.RR.txt2 = uicontrol(rrP1,'style','text',...
                  'string','Pausa Inicial [s]',...
                  'units','normalized','position',[.08 .487 .848 .094]);
              
%% Crear cuadro de Input

handles.RR.pausa = uicontrol(rrP1, 'style','edit',...
                     'units','normalized','position',[.08 .37 .848 .094]);
                                  
%% Crear mensajes de aviso

handles.RR.txt3 = uicontrol( rrP1, 'style','text',...
                      'string','Introduzca una pausa de hasta 10 segundos',...
                      'units','normalized','position',[.08 .04 .848 .138],...
                      'visible','on');               
                                                          
%% Crear botones de interacción

handles.RR.acep = uicontrol('parent',rrP1,...
                  'style','pushbutton',...
                  'string','Aceptar',...
                  'units','normalized','position',[.08 .214 .848 .094],...
                  'callback',{@RR_error,handles});

handles.RR.start = uicontrol('parent',rrP2,...
                  'style','togglebutton',...
                  'string','Comenzar/Detener',...
                  'units','normalized','position',[.08 .562 .848 .314],...
                  'callback',{@Start_RR,handles});

handles.RR.volver = uicontrol('parent',rrP2,...
                  'style','pushbutton',...
                  'string','Volver',...
                  'units','normalized','position',[.08 .162 .848 .314],...
                  'callback',@New_Inicio); 
              
end

%% VENTANA SINE SWEEP

function SSmedicion(~,~)
close all; 

%% Crear una figura para usar la GUI

f2 = figure('Visible','off','Position',[0 0 700 400],'Resize','off');
    movegui(f2,'center');
    f2.Name = 'Sine Sweep';
    f2.MenuBar = 'none';
    f2.NumberTitle = 'off';
    f2.Visible = 'on';  
        
%% Crear paneles

ssP1 = uipanel (f2,'title', '',...
                 'units','normalized','Position',[.04 .06 .6 .907]);
             
ssP2 = uipanel (f2,'title', '',...
                 'units','normalized','Position',[.65 .06 .31 .907]);
     
%% Crear título y títulos de Inputs

handles.SS.txt0 = uicontrol(ssP2,'style','text','FontSize',14,...
                  'string','Sine Sweep para medición del TR',...
                  'units','normalized','position',[.08 .63 .848 .183]);

handles.SS.txt1 = uicontrol(ssP1,'style','text',...
                       'string','Duración del Sine Sweep [s]',...
                       'units','normalized','position',[.073 .85 .380 .076]);
              
handles.SS.txt2 = uicontrol(ssP1,'style','text',...
                     'string','Pausa Inicial [s]',...
                     'units','normalized','position',[.550 .85 .380 .076]);
                
handles.SS.txt3 = uicontrol(ssP1,'style','text',...
                      'string','Frecuencia Inferior [Hz]',...
                      'units','normalized','position',[.076 .6 .380 .076]);
                
handles.SS.txt4 = uicontrol(ssP1,'style','text',...
                      'string','Frecuencia Superior [Hz]',...
                      'units','normalized','position',[.550 .6 .380 .076]);
                                   
handles.SS.txt5 = uicontrol(ssP1,'style','text',...
                         'string','Cantidad de Mediciones',...
                         'units','normalized','position',[0.076 0.35 0.380 0.076],...
                         'visible','on');                     

handles.SS.txt6 = uicontrol(ssP1,'style','text',...
                         'string','Cantidad de Posiciones',...
                         'units','normalized','position',[0.550 0.35 0.380 0.076],...
                         'visible','on');  
              
%% Crear mensajes de aviso y error

handles.SS.eD = uicontrol(ssP1,'style','text',...
                         'string','Entre 3 y 60 segundos',...
                         'units','normalized','position',[0.076 0.7 0.380 0.076],...
                         'visible','on');
                  
handles.SS.eP = uicontrol(ssP1,'style','text',...
                         'string','Entre 0 y 10 segundos',...
                         'units','normalized','position',[0.550 0.7 0.380 .076],...
                         'visible','on');
                     
handles.SS.eF = uicontrol(ssP1,'style','text',...
                         'string','Entre 20 y 20000 Hz',...
                         'units','normalized','position',[0.076 0.45 0.850 .076],...
                         'visible','on');

handles.SS.aviso = uicontrol(ssP1,'style','text',...
                         'string','Nota: de ser incompleto o inválido, el casillero se completará con un valor predeterminado',...
                         'units','normalized','position',[.07 .1 .88 .126],...
                         'visible','on');
            
%% Crear cuadros de Input

handles.SS.D = uicontrol(ssP1,'style','edit',...
                         'units','normalized','position',[.073 .8 .380 .076]);
 
handles.SS.P = uicontrol(ssP1,'style','edit',...
                       'string','0',...
                       'units','normalized','position',[.550 .8 .380 .076]);  
                     
handles.SS.F1 = uicontrol(ssP1,'style','edit',...
                        'units','normalized','position',[.073 .55 .380 .076]);
                
handles.SS.F2 = uicontrol(ssP1,'style','edit',...
                        'units','normalized','position',[.550 .55 .380 .076]);

handles.SS.R = uicontrol(ssP1,'style','popupmenu',...
                         'string',{'Cantidad de mediciones',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20},...
                         'max',1,'min',0,...
                         'units','normalized','position',[.073 .3 .380 .076]);
                
handles.SS.Pos = uicontrol(ssP1,'style','popupmenu',...
                         'string',{'Cantidad de posiciones',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20},...
                         'max',1,'min',0,...
                         'units','normalized','position',[.550 .3 .380 .076]);
                                    
%% Crear botones de interacción

handles.SS.acep = uicontrol(ssP1, 'style','pushbutton',...
                            'string','Aceptar',...
                            'units','normalized','position',[.650 .052 .280 .079],...
                            'callback',{@Error_SS,handles});

handles.SS.gen = uicontrol(ssP2,'style','pushbutton',...
                          'string','Generar Señal',...
                          'units','normalized','position',[.08 .45 .845 .094],...
                          'callback',{@Gen_SS,handles});

handles.SS.start = uicontrol(ssP2, 'style','pushbutton',...
                           'string','Comenzar',...
                           'units','normalized','position',[.08 .3 .845 .094],...
                           'callback',{@Start_SS,handles});             
                                  
handles.SS.volver = uicontrol(ssP2, 'style','pushbutton',...
                           'string','Volver',...
                           'units','normalized','position',[.08 .15 .845 .094],...
                           'callback',@New_Inicio); 

end

%% VENTANA CALCULO TR

function CalculoTR(~,~,handles)
close all; 

%% Crear una figura para usar la GUI
      
tr = figure('Visible','off','Position',[0 0 900 630],'Resize','off');
    movegui(tr,'north');
    tr.Name = 'Calculo TR';
    tr.MenuBar = 'none';
    tr.NumberTitle = 'off';
    tr.Visible = 'on';  
     
%% Crear paneles

handles.panel1 = uipanel (tr,'title', 'Input',...
                  'units','normalized','Position',[.008 .2 .197 .723]);  

handles.panel3 = uipanel (tr,'title', 'Grafico',...
                  'units','normalized','Position',[.204 .2 .590 .725]);

handles.panel4 = uipanel (tr,'title', 'Accion',...
                  'units','normalized','Position',[.793 .2 .197 .723]);    
              
handles.panel2 = uipanel (tr,'title', 'Descriptores promedio entre muestras',...
                  'units','normalized','Position',[.008 .016 .981 .191]);               

%% Titulos de Input   

handles.TR.txt1 = uicontrol(tr,'style','text','FontSize',14,...
                  'string','Cálculo del Tiempo de Reverberación',...
                  'units','normalized','position',[.02 .925 .96 .053]);

handles.TR.txt2 = uicontrol(handles.panel1,'style','text',...
                       'string','Seleccionar tipo de filtro',...
                       'units','normalized','position',[.097 .477 .826 .056]);  
                   
handles.TR.txt3 = uicontrol(handles.panel1,'style','text',...
                       'string','Seleccionar banda de freq que desee visualizar [Hz]',...
                       'units','normalized','position',[.097 .32 .826 .08]);                        
                   
handles.TR.txt4 = uicontrol(handles.panel1,'style','text',...
                       'string','Seleccionar medición que desee visualizar',...
                       'units','normalized','position',[.097 .16 .826 .08]);     
                   
handles.TR.txt5 = uicontrol(handles.panel4,'style','text','FontSize',10,...
                       'string','',...
                       'units','normalized','position',[.097 .69 .826 .08]);                    
                   
%% Preparar para calcular y calcular
                        
 handles.selectArchivo = uicontrol(handles.panel1, 'style','popupmenu',...
                            'string','Esperando archivos',...
                            'units','normalized','position',[.092 .107 .826 .056]);
   
 handles.table = uitable(handles.panel2,'units','normalized','Position',[.014 .078 .975 .878],...
                    'RowName',{'EDT','T10','T20','T30','C50','C80','D50','D80'},...
                    'ColumnName',{'b1','b2','b3','b4','b5','b6','b7','b8','b9','...'}); 
                
 handles.sinesweep = uicontrol(handles.panel1, 'style','pushbutton',...
                            'string','Generar Sine Sweep',...
                            'units','normalized','position',[.092 .72 .826 .107],...
                            'callback',{@SScalc,handles});
                        
 handles.abrir = uicontrol(handles.panel1, 'style','pushbutton',...
                            'string','Abrir archivos',...
                            'units','normalized','position',[.092 .853 .826 .107],...
                            'callback',{@abrir,handles});
                        
 handles.importarimpulso = uicontrol(handles.panel1, 'style','pushbutton',...
                            'string','Importar Impulso',...
                            'units','normalized','position',[.092 .587 .826 .107],...
                            'callback',{@abrirImp,handles});                       
                                             
 handles.selectBanda = uicontrol(handles.panel1, 'style','popupmenu',...
                            'string',{'31.5','63','125','250','500',...
                            '1000','2000','4000','8000','16000'},...
                            'units','normalized','position',[.092 .267 .826 .056]);                                                                                   
                                                               
 handles.selectFiltro = uicontrol(handles.panel1, 'style','popupmenu',...
                            'string',{'Por octava','Por tercio de octava'},...
                            'units','normalized','position',[.092 .424 .826 .056],...
                            'callback',{@selectFiltro,handles});
                        
%ventana de error
 handles.ventanaerror = figure('Position',[0 0 500 140],'Resize','off');
        movegui(handles.ventanaerror,'center');
        handles.ventanaerror.Name = 'ERROR';
        handles.ventanaerror.MenuBar = 'none';
        handles.ventanaerror.NumberTitle = 'off';
        handles.ventanaerror.Visible = 'off'; 
                        
        handles.error = uicontrol(handles.ventanaerror,'Style','text',...
            'Position',[50 20 400 100],...
            'String','Algo salió mal.',...
            'FontName','Trebuchet MS','FontSize',20); 
%                        
                       
 handles.calcular = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Calcular',...
                            'units','normalized','position',[.092 .82 .826 .14],...
                            'callback',{@Main,handles}); 
                       
 handles.resultados = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Ver resultados',...
                            'units','normalized','position',[.092 .562 .826 .091],...
                            'callback',{@table_Callback,handles}); 
                        
 handles.graficar = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Graficar',...
                            'units','normalized','position',[.092 .445 .826 .091],...
                            'callback',{@lookatthisgraph_Callback,handles}); 
                                                                  
 handles.guardar1 = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Guardar  Impulso',...
                            'units','normalized','position',[.092 .328 .826 .091],...
                            'callback',{@guardarimpulso,handles});
                        
 handles.guardar2 = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Guardar  resultados',...
                            'units','normalized','position',[.092 .211 .826 .091],...
                            'callback',{@guardartable,handles});
                       
 handles.volver = uicontrol(handles.panel4, 'style','pushbutton',...
                            'string','Volver',...
                            'units','normalized','position',[.092 .094 .826 .091],...
                            'callback', @New_Inicio);                           

%% crear grafico

handles.graf = axes(handles.panel3,'units','normalized','position', [.121 .129 .86 .835],...
                           'Xlim',[0 2],'Ylim', [-45,0]);
                    xlabel(handles.graf,'Tiempo [s]');
                    ylabel(handles.graf,'Amplitud [dBFS]');    
                    
       
end

%% FUNCIONES

% Funciones RR

function RR_error(~,~,handles)
% Funcion para el control de error de lo que se ingrese en la GUI 
% Si no se completa el casillero o el valor ingresado no es válido se
% corrige por un valor predeterminado.

global p
    p = str2double(get(handles.RR.pausa,'string'));
    handles.RR.txt3.Visible = 'Off';
    try 
      if p > 10 || isnan(p)
         error('Se espera un valor de tiempo de entre 0 y 10 segundos');
      end
    catch
         handles.RR.txt3.String = 'Se espera una pausa menor a 10 segundos';
         handles.RR.pausa.String = '5';
         handles.RR.txt3.Visible = 'On';
    end
end

function Start_RR(hObject,~,handles)
%   Funcion para reproducir y detener el ruido rosa.

    global p
    SoundCard = audioDeviceReader();
    handles.RR.txt3.Visible = 'Off';
    PinkNoise.SamplingFreq = SoundCard.SampleRate;
    PinkNoise.Duration=60;
    [x] = RuidoRosa(PinkNoise);
    PinkNoise.RR = audioplayer(x,PinkNoise.SamplingFreq);
    switch hObject.Value
        case 1
            handles.RR.txt3.String = sprintf('El ajuste comenzará en %d segundos.',p);
            handles.RR.txt3.Visible = 'On';
            pause(p);
            handles.RR.txt3.String = 'Reproduciendo ruido rosa para ajustar la fuente';
            play(PinkNoise.RR)
            handles.RR.txt3.Visible = 'On';
        case 0
            stop(PinkNoise.RR)
            handles.RR.txt3.Visible = 'Off';
    end
    setappdata(0,'PinkNoise',PinkNoise);
end

function [x]=RuidoRosa(PinkNoise)

    Nx = PinkNoise.Duration*PinkNoise.SamplingFreq;  % number of samples to synthesize
    B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
    A = [1 -2.494956002   2.017265875  -0.522189400];
    nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est. %roots raices de polinomio %round redondeo?
    v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
    x = filter(B,A,v);    % Apply 1/F roll-off to PSD
    x = x(nT60+1:end);    % Skip transient response
end

% Funciones SS

function Error_SS(~,~,handles)
% Funcion para el control de error de lo que se ingrese en la GUI 
% Si no se completan los casilleros o algun valor ingresado no es válido se
% corrigen por valores predeterminado.

%% Error del valor de la duracion del sine sweep
     handles.SS.eD.Visible = 'Off';
     SineSweep.Duration = str2double(get(handles.SS.D,'String'));
    try
       if SineSweep.Duration < 3 || SineSweep.Duration > 60 || isnan(SineSweep.Duration)
          error('Valor ingresado no válido.')
       end
    catch
        handles.SS.D.String = '5';
        SineSweep.Duration = str2double(get(handles.SS.D,'String'));
    end
    
%% Error del valor de la pausa
    handles.SS.eP.Visible = 'Off';
    SineSweep.Pause = str2double(get(handles.SS.P,'String'));
    try 
      if SineSweep.Pause > 10 || isnan(SineSweep.Pause)
         error('Valor ingresado no válido');
      end
    catch
         handles.SS.P = '5';
         SineSweep.Pause = str2double(get(handles.SS.P,'String'));
    end
    
 %% Errores de las frecuencias
    handles.SS.eF.Visible = 'Off';
    SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
    SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));
    try 
        if SineSweep.LowFreq<20||isnan(SineSweep.LowFreq)||SineSweep.LowFreq>=SineSweep.HighFreq||SineSweep.HighFreq>20000||isnan(SineSweep.HighFreq)
          error('Se esperan valores de frecuencia entre 20 y 20000 Hz. La frecuencia inferior no debe ser mayor que la superior.');
        end
    catch
        if SineSweep.LowFreq<20||isnan(SineSweep.LowFreq)
             handles.SS.F1.String = '250';
             SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
            
        end
        if SineSweep.HighFreq>20000||isnan(SineSweep.HighFreq)
            handles.SS.F2.String = '4000';
            SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));

        end
        if SineSweep.LowFreq>=SineSweep.HighFreq
            handles.SS.eF.String = 'La freq inferior no debe ser mayor que superior.'; 
            handles.SS.eF.Visible = 'On';
            handles.SS.F1.String = '250';
            SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
            handles.SS.F2.String = '4000';
            SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));
            
        end
    end
    
   %% Errores en # de medición
    handles.SS.avisoR.Visible = 'On';
    valmed = handles.SS.R.Value;
    medicion = handles.SS.R.String;
    SineSweep.med = medicion{valmed};
    try 
      if strcmp(SineSweep.med,'Cantidad de mediciones')
         error('La informacion ingresada no corresponde a un valor numérico');
      end
    catch
         handles.SS.R.Value = 6;
         valmed = handles.SS.R.Value;
         SineSweep.med = medicion{valmed};
    end
    
   %% Errores en # de posición
    handles.SS.avisoPos.Visible = 'on';
    valpos = handles.SS.Pos.Value;
    posicion = handles.SS.Pos.String;
    SineSweep.pos = posicion{valpos};
    try 
      if strcmp(SineSweep.pos,'Cantidad de posiciones')
         error('La informacion ingresada no corresponde a un valor numérico');
      end
    catch
         handles.SS.Pos.Value = 4 ;
         valpos = handles.SS.Pos.Value;
         SineSweep.pos = posicion{valpos};
    end
end

function Gen_SS(~,~,handles)
% Funcion para generar el sine sweep y el filtro inverso a partir de los
% valores ingresados (duración, frecuencia inicial, frecuencia final)
% así como también se guardan los valores de la cantidad de posiciones y
% mediciones en las que se llevará a cabo la práctica.

    handles.SS.aviso.Visible = 'off';
    SoundCard = audioDeviceReader();
    SineSweep.SamplingFreq = SoundCard.SampleRate;
    SineSweep.Duration = str2double(get(handles.SS.D,'String'));
    SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
    SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));
    SineSweep.pos = handles.SS.Pos.Value-1;
    SineSweep.med = handles.SS.R.Value-1;
    SineSweep.Pause = str2double(get(handles.SS.P,'String'));
    [SineSweep] = SineSweep_FiltroInverso(SineSweep);
    handles.SS.aviso.String = 'El sine sweep y el filtro inverso han sido generados.';
    handles.SS.aviso.Visible = 'on';
    setappdata(0,'SineSweep',SineSweep);
end

function Start_SS(~,~,handles)
% Funcion para realizar la medición del TR con el Sine Sweep.
% Realiza el total de mediciones pedidas en cada posición de manera
% continua, luego solicita informar el cambio de posición en caso de
% elegir más de una para continuar.
      
    global i j
    SineSweep = getappdata(0,'SineSweep');
    handles.SS.aviso.Visible = 'off';
    
    %mensajes de aviso que seran utilizados durante la medicion
    AP = sprintf('La medición comenzará en %d segundos.',SineSweep.Pause); 
    AR = sprintf('Reproduciendo sine sweep y obteniendo respuesta de la sala...');
    AN = sprintf('Mediciones en la posición %d completadas. Cambie de posición para continuar',i);
   
    i = 1; %posicion
    j = 1; %n° medicion en posicion i
    while i <= SineSweep.pos %medicion del TR
        while j <= SineSweep.med
            handles.SS.aviso.String = AP;
            handles.SS.aviso.Visible = 'on';
            pause(SineSweep.Pause);
            handles.SS.aviso.Visible = 'off';
            handles.SS.aviso.String = AR;
            handles.SS.aviso.Visible = 'on';
            [y] = Rep_Grab(SineSweep,i,j);
            j=j+1;
        end

        handles.SS.aviso.String = AN;
        handles.SS.aviso.Visible = 'on';
        j=1;
        resp = 0;
        if i < SineSweep.med
            handles.SS.aviso.Visible = 'off';
            AM = sprintf('Cambio de posicion registrado. Posición número %d. Para continuar con las mediciones presione Enter en Command Window.',i); %aviso de cambio de posicion
            handles.SS.aviso.String = AM;
            handles.SS.aviso.Visible = 'on';
            pause;
            i = i+1;
        else
            handles.SS.aviso.Visible = 'off';
            AM = sprintf('Última posición. Posición número %d',i);
            handles.SS.aviso.String = AM;
            handles.SS.aviso.Visible = 'on';
            break
        end
    end
end

function [resp] = Pos_SS

    resp = 1;
    
end

function [SineSweep]=SineSweep_FiltroInverso(SineSweep)

    w1 = 2*pi*SineSweep.LowFreq; %frecuencia angular inferior
    w2 = 2*pi*SineSweep.HighFreq; %frecuencia angular superior
    K = (SineSweep.Duration*w1)/log(w2/w1);
    t = 0:1/SineSweep.SamplingFreq:SineSweep.Duration;
    L = SineSweep.Duration/log(w2/w1);
    SS = sin(K*(exp(t/L)-1)); %esto es el sinesweep
    SineSweep.Sweep = SS/abs(max(SS));
    %filtro inverso
    w = (K/L)*exp(t/L);
    m=w1./(2*pi*w); %modulacion
    FI = m.*(fliplr(SineSweep.Sweep)); %filtro inverso
    SineSweep.Filter = FI/abs(max(FI)); %normalizacion amplitud
    SineSweep.Interval = t;
end

function [y] = Rep_Grab(SineSweep,i,j)

    SSRec=audiorecorder(SineSweep.SamplingFreq,16,1);
    SS1=audioplayer(SineSweep.Sweep,SineSweep.SamplingFreq);
    tic;
    tin=tic;
    record(SSRec,SineSweep.Duration+2);
    play(SS1);
    tfin=toc(tin);
    sound(SineSweep.Sweep,SineSweep.SamplingFreq);
    pause(SineSweep.Duration+3);
    y = getaudiodata(SSRec);
    %disp(tfin);
    filename = sprintf('Posicion%d_Medicion%d.wav',i,j);
    audiowrite(filename,y,SineSweep.SamplingFreq);
end

% Funciones Cálculo TR

function abrir(~,~,handles)
% Funcion para abrir los archivos de las mediciones realizadas. 
% permite elegir multiples archivos a la vez
% permite abrir archivos en formato .wav o .mat 
    
    [filename, pathname] = uigetfile({'*.wav';'*.mat'},...
        'Seleccione todas las mediciones','Multiselect','on');
    pathfile = strcat(pathname,filename);
    a = ischar(filename);
    if a == 1
        set(handles.selectArchivo,'String',filename);
        [matrix.medicion.m1,matrix.SamplingFreq] = audioread(pathfile);
        matrix.n = 1;
    else
        set(handles.selectArchivo,'String',char(filename));
        [x,fs] = cellfun(@audioread,pathfile,'UniformOutput',0);
        matrix.n = numel(x);
        for i = 1:matrix.n
            matrix.medicion.(sprintf('m%d',i)) = x{i};
        end
    matrix.SamplingFreq = fs{1};
    end
    matrix.op = 0;
    setappdata(0,'matrix',matrix);
end

function abrirImp(~,~,handles)
% Funcion para abrir archivos de impulsos, en lugar de archivos de
% mediciones. Para realizar el análsis de los descriptores del mismo.
% permite elegir multiples archivos a la vez. 
% permite abrir archivos en formato .wav o .mat 
    
    [filename, pathname] = uigetfile({'*.wav';'*.mat'},...
        'Seleccione todas las mediciones','Multiselect','on');
    pathfile = strcat(pathname,filename);
    a = ischar(filename);
    if a == 1
        set(handles.selectArchivo,'String',filename);
        [matrix.RI.m1,matrix.SamplingFreq] = audioread(pathfile);
        matrix.n = 1;
    else
        set(handles.selectArchivo,'String',char(filename));
        [x,fs] = cellfun(@audioread,pathfile,'UniformOutput',0);
        matrix.n = numel(x);
        for i = 1:matrix.n
            matrix.RI.(sprintf('m%d',i)) = x{i};
        end
    matrix.SamplingFreq = fs{1};
    end
    matrix.op = 1;
    setappdata(0,'matrix',matrix);
end

function SScalc(~,~,handles)
% SScalc permite obtención de parametros necesarios 
% (frecuencia inicial, frecuencia final y duración) 
% para generar el barrido con el que se realizará la convolución 
% si no se ha generado un Sine Sweep en la ventana de medición 

    SSParam = figure('Visible','on','Position',[0 0 400 250],'Resize','off');
    movegui(SSParam,'center');
    SSParam.Name = 'Sine Sweep para medición del TR';
    SSParam.MenuBar = 'none';
    SSParam.NumberTitle = 'off';
    
    ssP1 = uipanel (SSParam,'title', 'Input',...
                 'units','normalized','Position',[.055 .060 .9 .90]);
    
    handles.SScalc.txt1 = uicontrol(ssP1,'style','text',...
                       'string','Duración del Sine Sweep [s]',...
                       'units','normalized','position',[.3 .896 .41 .08]);
    handles.SScalc.txt2 = uicontrol(ssP1,'style','text',...
                      'string','Frecuencia Inferior [Hz]',...
                      'units','normalized','position',[.076 .570 .380 .08]);             
    handles.SScalc.txt3 = uicontrol(ssP1,'style','text',...
                      'string','Frecuencia Superior [Hz]',...
                      'units','normalized','position',[.550 .570 .380 .08]);
     
    handles.SScalc.txt4 = uicontrol(ssP1,'style','text',...
                      'string','Nota: de ser incompleto o inválido, el casillero se completará con un valor predeterminado',...
                      'units','normalized','position',[.076 .18 .85 .18]);
                  
    handles.SS.eD = uicontrol(ssP1,'style','text',...
                         'string','Entre 3 y 60 segundos',...
                         'units','normalized','position',[0.3 0.700 0.41 0.08],...
                         'visible','on');

    handles.SS.eF = uicontrol(ssP1,'style','text',...
                         'string','Entre 20 y 20000 Hz',...
                         'units','normalized','position',[0.076 0.370 0.850 .08],...
                         'visible','on');                  
                     
    handles.SS.D = uicontrol(ssP1,'style','edit',...
                          'units','normalized','position',[.3 .800 .41 .08]);

    handles.SS.F1 = uicontrol(ssP1,'style','edit',...
                          'units','normalized','position',[.073 .480 .380 .08]);

    handles.SS.F2 = uicontrol(ssP1,'style','edit',...
                          'units','normalized','position',[.550 .480 .380 .08]);

    handles.SS.gen = uicontrol(ssP1,'style','pushbutton',...
                          'string','Generar Señal',...
                          'units','normalized','position',[.073 .07 .380 .15],...
                          'callback',{@Gen_calcSS,handles});   

    handles.SS.close = uicontrol(ssP1,'style','pushbutton',...
                          'string','Volver',...
                          'units','normalized','position',[.550 .07 .380 .15],...
                          'callback', 'close(gcf)');                       
end

function selectFiltro(object_handle,~,handles)
% selectFiltro determina las opciones de frecuencias centrales de los anchos de
% banda posibles, a partir de la selección de filtro por octava o tercio
% de octava, para la posterior visualización de los resultados obtenidos.

    try
        switch get(object_handle,'Value')
            case 1
            set(handles.selectBanda,'Value',1);
            set(handles.selectBanda,'String',{'31.5','63','125',...
                '250','500','1000','2000','4000','8000','16000'});
            case 2
            set(handles.selectBanda,'Value',1);
            set(handles.selectBanda,'String',{'28.5','32','36','57',...
                '64','71.8','111.3','125','140.3','222.7','250',...
                '280.6','445.5','500','561','890.9','1000','1122.46',...
                '1781.79','2000','2244.92','3563.59','4000','4489.84',...
                '7127.19','8000','8979.7','14254.4','16000','17959.4'});
        end
    catch
            set(handles.selectBanda,'Value',1);
            set(handles.selectBanda,'String',{'31.5','63','125',...
                '250','500','1000','2000','4000','8000','16000'});
    end
end

function Main(~,~,handles)
% Main realiza los todos los cálculos necesarios para la obtención de los
% descriptores y el tiempo de reverberación

    handles.TR.txt5.Visible = 'Off';
     try
         if strcmp(handles.selectArchivo,'Esperando archivos')
             error('No seleccionaste nada pa');
         end
         
    SineSweep = getappdata(0,'SineSweep');
    matrix = getappdata(0,'matrix');
    if matrix.op == 0
        matrix = RespuestaImpulsiva(matrix,SineSweep);
    elseif matrix.op == 1
        matrix.RI = structfun(@T_Inicial,matrix.RI,'UniformOutput',0);
    end
    [T] = lundeby(matrix);
    for i = 1:matrix.n
        matrix.RI.(sprintf('m%d',i)) = matrix.RI.(sprintf('m%d',i))(1:T);
    end
    
    handles.TR.txt5.String = 'Calculando...';
    handles.TR.txt5.Visible = 'On';
    matrix = Filtros(matrix);
    matrix = Plotter(matrix);
    matrix = Suaaaveee(matrix);
    matrix = Ct0(matrix);
    matrix = Dt0(matrix);
    matrix = dBfs(matrix);
    matrix = TR(matrix);
    [v] = Convert(matrix);
    [v] = Promedio(v);
    [v] = Tablatize(v);
    setappdata(0,'matrix',matrix);
    setappdata(0,'v',v);
    handles.TR.txt5.String = 'Descriptores calculados exitosamente :)';
     catch
        handles.ventanaerror.Visible = 'on';
        pause(1)
        set(handles.ventanaerror,'Visible','off');    
     end
end

function lookatthisgraph_Callback(~,~,handles)
% Función para graficar los resultados obtenidos a partir de la selección
% de un archivo de medición y una banda de frecuencia especificas.

    try
    matrix = getappdata(0,'matrix');
    v = getappdata(0,'v');
    i = get(handles.selectFiltro,'Value');
    j = get(handles.selectBanda,'Value');
    k = get(handles.selectArchivo,'Value');

    if i==2
        i=3;
    end

    fs = matrix.SamplingFreq;
    x    = matrix.(sprintf('plotter%d',i)).(sprintf('b%dm%d',j,k));
    sch  = matrix.(sprintf('dB%d',i)).(sprintf('b%dm%d',j,k));
    n    = 1/fs:1/fs:length(x)/fs;
    nEDT = matrix.(sprintf('nEDT%d',i)).(sprintf('b%dm%d',j,k));
    pEDT = matrix.(sprintf('pEDT%d',i)).(sprintf('b%dm%d',j,k));
    rEDT = pEDT(1)*nEDT+pEDT(2);
    nT10 = matrix.(sprintf('nT10%d',i)).(sprintf('b%dm%d',j,k));
    pT10 = matrix.(sprintf('pT10%d',i)).(sprintf('b%dm%d',j,k));
    rT10 = pT10(1)*nT10+pT10(2);
    nT20 = matrix.(sprintf('nT20%d',i)).(sprintf('b%dm%d',j,k));
    pT20 = matrix.(sprintf('pT20%d',i)).(sprintf('b%dm%d',j,k));
    rT20 = pT20(1)*nT20+pT20(2);
    nT30 = matrix.(sprintf('nT30%d',i)).(sprintf('b%dm%d',j,k));
    pT30 = matrix.(sprintf('pT30%d',i)).(sprintf('b%dm%d',j,k));
    rT30 = pT30(1)*nT30+pT30(2);

    cla %para limpiar ejes

    plot(n,x,':c','LineWidth',1);
    hold on
    plot(n,sch,'--r','LineWidth',1);
    xlim([n(1),nT30(end)]);
    ylim([-45,0]);
    plot(nEDT,rEDT,'-.b','LineWidth',1);
    plot(nT10,rT10,'-.k','LineWidth',1);
    plot(nT20,rT20,'--m','LineWidth',1);
    plot(nT30,rT30,':g','LineWidth',1);

    legend('Envolvente','Schroeder Decay','EDT','T10','T20','T30');
    xlabel('Tiempo [s]');
    ylabel('Amplitud [dBFS]');
    catch
        handles.ventanaerror.Visible = 'on';
        pause(1)
        set(handles.ventanaerror,'Visible','off');
    end
    
end

function table_Callback(~,~,handles)
% Función para devolver los resultados de los descriptores, en cada
% frecuencia central, según se seleccione el filtro por octavas o tercios
% de octavas. 
% Si se abre un sólo archivo entrega los resultados
% correspondientes a esa medición
% Si se abren múltiples archivos entrega un promedio de todas las
% mediciones.

    try
    matrix = getappdata(0,'matrix');
    v = getappdata(0,'v');
    i = get(handles.selectFiltro,'Value');
    j = get(handles.selectBanda,'Value');
    k = get(handles.selectArchivo,'Value');

    if i==2
        i=3;
    end

    fs = matrix.SamplingFreq;
    x    = matrix.(sprintf('plotter%d',i)).(sprintf('b%dm%d',j,k));
    sch  = matrix.(sprintf('dB%d',i)).(sprintf('b%dm%d',j,k));
    n    = 1/fs:1/fs:length(x)/fs;
    nEDT = matrix.(sprintf('nEDT%d',i)).(sprintf('b%dm%d',j,k));
    pEDT = matrix.(sprintf('pEDT%d',i)).(sprintf('b%dm%d',j,k));
    rEDT = pEDT(1)*nEDT+pEDT(2);
    nT10 = matrix.(sprintf('nT10%d',i)).(sprintf('b%dm%d',j,k));
    pT10 = matrix.(sprintf('pT10%d',i)).(sprintf('b%dm%d',j,k));
    rT10 = pT10(1)*nT10+pT10(2);
    nT20 = matrix.(sprintf('nT20%d',i)).(sprintf('b%dm%d',j,k));
    pT20 = matrix.(sprintf('pT20%d',i)).(sprintf('b%dm%d',j,k));
    rT20 = pT20(1)*nT20+pT20(2);
    nT30 = matrix.(sprintf('nT30%d',i)).(sprintf('b%dm%d',j,k));
    pT30 = matrix.(sprintf('pT30%d',i)).(sprintf('b%dm%d',j,k));
    rT30 = pT30(1)*nT30+pT30(2);

    switch i
        case 1
            set(handles.table,'ColumnName',{'31.5','63','125','250',...
                '500','1000','2000','4000','8000','16000'});
        case 3
            set(handles.table,'ColumnName',{'28.5','32','36','57','64',...
                '71.8','111.3','125','140.3','222.7','250','280.6','445.5',...
                '500','561','890.9','1000','1122.46','1781.79','2000',...
                '2244.92','3563.59','4000','4489.84','7127.19','8000',...
                '8979.7','14254.4','16000','17959.4'});
    end
    set(handles.table,'Data',v.(sprintf('Tabla%d',i)));
    catch
        handles.ventanaerror.Visible = 'on';
        pause(1)
        set(handles.ventanaerror,'Visible','off');
    end
end

function [matrix] = RespuestaImpulsiva(matrix,SineSweep)

    [matrix.RI] = structfun(@RI,matrix.medicion,'UniformOutput',0);
    function [h] = RI(y)
        y = transpose(y);
        [SineSweep] = SineSweep_FiltroInverso(SineSweep);
        n1 = numel(SineSweep.Filter);
        n2 = numel(y);
        SineSweep.Filter = [SineSweep.Filter zeros(1,n1-1)];
        y = [y zeros(1,n2-1)];
        l1 = length(y);
        l2 = length(SineSweep.Filter);
        SineSweep.Filter = [SineSweep.Filter zeros(1,l1-l2)];
        fy = fft(y);
        fk = fft(SineSweep.Filter);
        fh = fy.*fk;
        h = ifft(fh);
        h = h/max(abs(h));

        [h] = T_Inicial(h);
    end
end

function [h] = T_Inicial(h)

    maximo = find(abs(h) == max(abs(h)));
    energia = (h/h(maximo(1))).^2;
    aux = energia(1:maximo(1)-1);
    inicio = maximo(1);
    if any(aux > 0.01)
        inicio = find(aux > 0.01, 1 );
        aux = energia(1:inicio-1);
    end
    h = h(inicio:end);
end

function [truncamiento] = lundeby(matrix)

    h = matrix.RI.m1;
    energia = h.^2;

    rms_dB = 10*log10(mean(energia(round(.9*length(energia)):end))/max(energia));
    t = floor(length(energia)/matrix.SamplingFreq/0.32);
    v = floor(length(energia)/t);
    for n=1:t
        media(n) = mean(energia((((n-1)*v)+1):(n*v)));
        tiempo(n) = ceil(v/2)+((n-1)*v);
    end
    mediadB = 10*log10(media/max(energia));

    r = find(mediadB > rms_dB+10, 1, 'last' );
    if any (mediadB(1:r) < rms_dB+10)
        r = find(mediadB(1:r) < rms_dB+10, 1 );
    end
    if isempty(r)
        r=10;
    elseif r<10
        r=10;
    end
    [A,B] = RegresionLineal(tiempo(1:r),mediadB(1:r));
    CR = (rms_dB-A)/B;
    if rms_dB > -20
        truncamiento=length(energia);
    else
        error=1;
        INTMAX=50;
        veces=1;
        while (error > 0.0001 && veces <= INTMAX)
            clear r t v n media tiempo;
            p = 5;
            delta = abs(10/B);
            v = floor(delta/p);
            t = floor(length(energia(1:round(CR-delta)))/v);
            if t < 2
                t=2;
            elseif isempty(t)
                t=2;
            end

            for n=1:t
                media(n) = mean(energia((((n-1)*v)+1):(n*v)));
                tiempo(n) = ceil(v/2)+((n-1)*v);
            end
            mediadB = 10*log10(media/max(energia));

            clear A B noise rms_dB;
            [A,B] = RegresionLineal(tiempo,mediadB);

            noise = energia(round(CR+delta):end);
            if (length(noise) < round(.1*length(energia)))
                noise = energia(round(.9*length(energia)):end); 
            end       
            rms_dB = 10*log10(mean(noise)/max(energia));

            error = abs(CR - (rms_dB-A)/B)/CR;
            CR = round((rms_dB-A)/B);
            veces = veces + 1;
        end
    end

        if CR > length(energia)
            truncamiento = length(energia);
        else
            truncamiento = CR;
        end
end

function [matrix] = Filtros(matrix)

    h = matrix.RI;
    n = length(fieldnames(h));
    for j = 1:n
         x = h.(sprintf('m%d',j));
        Oct = [31.62 63.1 125.89 251.19 501.19 1000 1995.26 3981.07 7943.28 15848.93];
            for i = 1:length(Oct)
                filter1 = fdesign.octave(1,'Class 0','N,F0',6,Oct(i),matrix.SamplingFreq);
                filter1 = design(filter1,'Butter');
                matrix.Oct.(sprintf('b%dm%d',i,j)) = filter(filter1,x);
            end

        Ter_Oct =[28.5 32 36 57 64 71.8 111.3 125 140.3 222.7 250 280.6 445.5 500 561 890.9 1000 1122.46 1781.79 2000 2244.92 3563.59 4000 4489.84 7127.19 8000 8979.7 14254.4 16000 17959.4];
            for i = 1:length(Ter_Oct)
                filter1 = fdesign.octave(3,'Class 0','N,F0',8,Ter_Oct(i),matrix.SamplingFreq);
                filter1 = design(filter1,'Butter');
                matrix.Ter_Oct.(sprintf('b%dm%d',i,j)) = filter(filter1,x);
            end  
    end
end

function [matrix] = Suaaaveee(matrix)

matrix.SuaveOct = structfun(@Suavizado,matrix.Oct,'UniformOutput',0);
matrix.SuaveTer_Oct = structfun(@Suavizado,matrix.Ter_Oct,'UniformOutput',0);

    function y = Suavizado(x)
    a = abs(hilbert(x));
    y = cumsum(a.^2,'reverse');
    end
end

function [matrix] = Ct0(matrix)

    matrix.C50Oct = structfun(@calcC50,matrix.SuaveOct,'UniformOutput',0);
    matrix.C50Ter_Oct = structfun(@calcC50,matrix.SuaveTer_Oct,'UniformOutput',0);
    
    matrix.C80Oct = structfun(@calcC80,matrix.SuaveOct,'UniformOutput',0);
    matrix.C80Ter_Oct = structfun(@calcC80,matrix.SuaveTer_Oct,'UniformOutput',0);
    
    function y = calcC50(x)
        SoundCard = audioDeviceReader();
        fs = SoundCard.SampleRate;
        n = length(x);
        y1 = x(1:0.05*fs).^2;
        y2 = x(0.05*fs:n).^2;
        y = 10*log10(sum(y1)/sum(y2));
    end

    function y = calcC80(x)
        SoundCard = audioDeviceReader();
        fs = SoundCard.SampleRate;
        n = length(x);
        y1 = x(1:0.08*fs).^2;
        y2 = x(0.08*fs:n).^2;
        y = 10*log10(sum(y1)/sum(y2));
    end
end

function [matrix] = Dt0(matrix)

    matrix.D80Oct = structfun(@calcD80,matrix.SuaveOct,'UniformOutput',0);
    matrix.D80Ter_Oct = structfun(@calcD80,matrix.SuaveTer_Oct,'UniformOutput',0);
    
    matrix.D50Oct = structfun(@calcD50,matrix.SuaveOct,'UniformOutput',0);
    matrix.D50Ter_Oct = structfun(@calcD50,matrix.SuaveTer_Oct,'UniformOutput',0);
    
    function y = calcD80(x)
        SoundCard = audioDeviceReader();
        fs = SoundCard.SampleRate;
        n = length(x);
        y1 = x(1:0.08*fs).^2;
        y2 = x(1:n).^2;
        y = 100*(sum(y1)/sum(y2));
    end

    function y = calcD50(x)
        SoundCard = audioDeviceReader();
        fs = SoundCard.SampleRate;
        n = length(x);
        y1 = x(1:0.05*fs).^2;
        y2 = x(1:n).^2;
        y = 100*(sum(y1)/sum(y2));
    end
end

function [matrix] = dBfs(matrix)

matrix.dB1 = structfun(@dBFull_Scale,matrix.SuaveOct,'UniformOutput',0);
matrix.dB3 = structfun(@dBFull_Scale,matrix.SuaveTer_Oct,'UniformOutput',0);

    function y = dBFull_Scale(x)
        y = 10*log10(x/max(x));
    end
end

function matrix = TR(matrix)

fs = matrix.SamplingFreq;

[matrix.nEDT1, matrix.pEDT1, matrix.edt1, matrix.nT101, matrix.pT101, matrix.t101,...
    matrix.nT201, matrix.pT201, matrix.t201, matrix.nT301, matrix.pT301, matrix.t301] =...
    structfun(@TRc,matrix.dB1,'UniformOutput',0);
[matrix.nEDT3, matrix.pEDT3, matrix.edt3, matrix.nT103, matrix.pT103, matrix.t103,...
    matrix.nT203, matrix.pT203, matrix.t203, matrix.nT303, matrix.pT303, matrix.t303] =...
    structfun(@TRc,matrix.dB3,'UniformOutput',0);

    function [nEDT,pEDT,edt,nT10,pT10,t10,nT20,pT20,t20,nT30,pT30,t30] = TRc(x)

        noise = find(x<=-45,1,'first');     
        x = x(1:noise);
        n = 1/fs:1/fs:length(x)/fs;                 %Dominio Temporal
        
        %EDT
        sup = find(x==0,1,'first');
        low = find(x<=-10,1,'first');
        nEDT = n(sup:low);
        xEDT = x(sup:low);
        pEDT = RegresionLineal_Matriz(nEDT,xEDT);
        edt = (-60-pEDT(2))/pEDT(1);
        
        %T10
        sup = find(x<=-5,1,'first');   
        low = find(x<=-15,1,'first');
        nT10 = n(sup:low);
        xT10 = x(sup:low);
        pT10 = RegresionLineal_Matriz(nT10,xT10);
        t10 = (-60-pT10(2))/pT10(1);
        
        %T20
        sup = find(x<=-5,1,'first');   
        low = find(x<=-25,1,'first');
        nT20 = n(sup:low);
        xT20 = x(sup:low);
        pT20 = RegresionLineal_Matriz(nT20,xT20);
        t20 = (-60-pT20(2))/pT20(1);
        
        %T30
        sup = find(x<=-5,1,'first');   
        low = find(x<=-35,1,'first');
        nT30 = n(sup:low);
        xT30 = x(sup:low);
        pT30 = RegresionLineal_Matriz(nT30,xT30);
        t30 = (-60-pT30(2))/pT30(1);
        
    end
end

function v = Convert(matrix)

    if matrix.op == 0
        n = length(fieldnames(matrix.medicion));
    elseif matrix.op == 1
        n = length(fieldnames(matrix.RI));
    end

    for i = 1:n
        for j = 1:10
            v.edt1.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.edt1.(sprintf('b%dm%d',j,i));
            v.t101.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t101.(sprintf('b%dm%d',j,i));
            v.t201.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t201.(sprintf('b%dm%d',j,i));
            v.t301.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t301.(sprintf('b%dm%d',j,i));
            v.C501.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.C50Oct.(sprintf('b%dm%d',j,i));
            v.C801.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.C80Oct.(sprintf('b%dm%d',j,i));
            v.D501.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.D50Oct.(sprintf('b%dm%d',j,i));
            v.D801.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.D80Oct.(sprintf('b%dm%d',j,i));
        end
        for j = 1:30
            v.edt3.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.edt3.(sprintf('b%dm%d',j,i));
            v.t103.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t103.(sprintf('b%dm%d',j,i));
            v.t203.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t203.(sprintf('b%dm%d',j,i));
            v.t303.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.t303.(sprintf('b%dm%d',j,i));
            v.C503.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.C50Ter_Oct.(sprintf('b%dm%d',j,i));
            v.C803.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.C80Ter_Oct.(sprintf('b%dm%d',j,i));
            v.D503.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.D50Ter_Oct.(sprintf('b%dm%d',j,i));
            v.D803.(sprintf('b%d',j)).(sprintf('m%d',i)) = matrix.D80Ter_Oct.(sprintf('b%dm%d',j,i));
        end
    end
end

function [v] = Promedio(v)

    for i = 1:10
    V = struct2cell(v.edt1.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.edt1.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t101.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t101.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t201.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t201.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t301.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t301.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.C501.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.C501.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.C801.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.C801.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.D501.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.D501.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.D801.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.D801.(sprintf('b%d',i)).avg = mean(C);
    end

    for i = 1:30
    V = struct2cell(v.edt3.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.edt3.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t103.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t103.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t203.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t203.(sprintf('b%d',i)).avg = mean(C);

    V = struct2cell(v.t303.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.t303.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.C503.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.C503.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.C803.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.C803.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.D503.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.D503.(sprintf('b%d',i)).avg = mean(C);
    
    V = struct2cell(v.D803.(sprintf('b%d',i)));
    C = cat(1,V{:});
    v.D803.(sprintf('b%d',i)).avg = mean(C);
    end
end

function [v] = Tablatize(v)

    edt1 = zeros(10,1);
    t101 = zeros(10,1);
    t201 = zeros(10,1);
    t301 = zeros(10,1);
    c501 = zeros(10,1);
    c801 = zeros(10,1);
    d501 = zeros(10,1);
    d801 = zeros(10,1);
    edt3 = zeros(30,1);
    t103 = zeros(30,1);
    t203 = zeros(30,1);
    t303 = zeros(30,1);
    c503 = zeros(30,1);
    c803 = zeros(30,1);
    d503 = zeros(30,1);
    d803 = zeros(30,1);

    for i = 1:10
        edt1(i)=v.edt1.(sprintf('b%d',i)).avg;
        t101(i)=v.t101.(sprintf('b%d',i)).avg;
        t201(i)=v.t201.(sprintf('b%d',i)).avg;
        t301(i)=v.t301.(sprintf('b%d',i)).avg;
        c501(i)=v.C501.(sprintf('b%d',i)).avg;
        c801(i)=v.C801.(sprintf('b%d',i)).avg;
        d501(i)=v.D501.(sprintf('b%d',i)).avg;
        d801(i)=v.D801.(sprintf('b%d',i)).avg;
    end

    for i = 1:30
        edt3(i)=v.edt3.(sprintf('b%d',i)).avg;
        t103(i)=v.t103.(sprintf('b%d',i)).avg;
        t203(i)=v.t203.(sprintf('b%d',i)).avg;
        t303(i)=v.t303.(sprintf('b%d',i)).avg;
        c503(i)=v.C503.(sprintf('b%d',i)).avg;
        c803(i)=v.C803.(sprintf('b%d',i)).avg;
        d503(i)=v.D503.(sprintf('b%d',i)).avg;
        d803(i)=v.D803.(sprintf('b%d',i)).avg;
    end

    TR1 = cat(2,edt1,t101,t201,t301,c501,c801,d501,d801);
    TR3 = cat(2,edt3,t103,t203,t303,c503,c803,d503,d803);

    v.Tabla1 = TR1';
    v.Tabla3 = TR3';

    fields = {'EDT';'T10';'T20';'T30';'C50';'C80';'D50';'D80'};
    
    f31_5 = v.Tabla1(:,[1]); f63 = v.Tabla1(:,[2]); f125 = v.Tabla1(:,[3]);
    f250 = v.Tabla1(:,[4]); f500 = v.Tabla1(:,[5]); f1000 = v.Tabla1(:,[6]);
    f2000 = v.Tabla1(:,[7]); f4000 = v.Tabla1(:,[8]); f8000 = v.Tabla1(:,[9]);
    f16000 = v.Tabla1(:,[10]);
    
    v.T1 = table(f31_5,f63,f125,f250,f500,f1000,f2000,f4000,f8000,f16000,...
                                    'RowNames',fields);
                                
    f28_5 = v.Tabla3(:,[1]); f32 = v.Tabla3(:,[2]); f36 = v.Tabla3(:,[3]);
    f57 = v.Tabla3(:,[4]); f64 = v.Tabla3(:,[5]); f71_8 = v.Tabla3(:,[6]);
    f111_3 = v.Tabla3(:,[7]); f125 = v.Tabla3(:,[8]); f140_3 = v.Tabla3(:,[9]);
    f222_7 = v.Tabla3(:,[10]); f250 = v.Tabla3(:,[11]); f280_6 = v.Tabla3(:,[12]);
    f445_5 = v.Tabla3(:,[13]); f500 = v.Tabla3(:,[14]); f561 = v.Tabla3(:,[15]);
    f890_9 = v.Tabla3(:,[16]); f1000 = v.Tabla3(:,[17]); f1122_46 = v.Tabla3(:,[18]);
    f1781_79 = v.Tabla3(:,[19]); f2000 = v.Tabla3(:,[20]); f2244_92 = v.Tabla3(:,[21]);
    f3563_59 = v.Tabla3(:,[22]); f4000 = v.Tabla3(:,[23]); f4489_84 = v.Tabla3(:,[24]);
    f7127_19 = v.Tabla3(:,[25]); f8000 = v.Tabla3(:,[26]); f8979 = v.Tabla3(:,[27]);
    f14254_4 = v.Tabla3(:,[28]); f16000 = v.Tabla3(:,[29]); f17959_4 = v.Tabla3(:,[30]);
    
    v.T3 = table(f28_5,f32,f36,f57,f64,f71_8,f111_3,f125,f140_3,f222_7,f250,...
        f280_6,f445_5,f500,f561,f890_9,f1000,f1122_46,f1781_79,f2000,f2244_92,...
        f3563_59,f4000,f4489_84,f7127_19,f8000,f8979,f14254_4,f16000,f17959_4,...
                                                        'RowNames',fields);
end

function [a0,a1] = RegresionLineal (x,y)

    a0 = mean(y) - ((sum(x.*y)-(mean(y)*sum(x)))/(sum(x.^2)-(mean(x)*sum(x))))*mean(x);
    a1 = (sum(x.*y) -(mean(y)*(sum(x))))/(sum(x.^2) - (mean(x)*sum(x)));

end

function [p] = RegresionLineal_Matriz(n,x)       

    N=length(n); %numero de elementos que contiene n
        I=ones(N,1);
        Z=[n(:) I(:)];
        c=Z'*Z;
        o=Z'*(x(:));
        p=c\o; %p(1)=pendiente p(2)=oo
        promedion=sum(n)/N;
        promediox=sum(x)/N;
        %calculo del error r=coeficiente de correlacion
%       Snx=(sum(n.*x)/N)-(promedion*promediox);
%       Sn=sqrt((sum(n.^2)/N)-(sum(n)/N)^2);
%       Sx=sqrt((sum(x.^2)/N)-(sum(x)/N)^2);
%       r=(Snx/(Sn*Sx));
        %p(3) = r^2;   %porcentaje de exactitud de la regresion lineal
end

function matrix = Plotter(matrix)

    matrix.plotter1 = structfun(@plotit,matrix.Oct,'UniformOutput',0);
    matrix.plotter3 = structfun(@plotit,matrix.Ter_Oct,'UniformOutput',0);

    function y = plotit(x)
        a = abs(hilbert(x));
        y = 20*log10(a/max(a));
    end
end

function Gen_calcSS(~,~,handles)
% funcion para generar el sine sweep y el filtro inverso  

    SoundCard = audioDeviceReader();
    SineSweep.SamplingFreq = SoundCard.SampleRate;
    handles.SScalc.txt4.Visible = 'off';

%% Error de la duración
    SineSweep.Duration = str2double(get(handles.SS.D,'String'));    
    try
       if SineSweep.Duration < 3 || SineSweep.Duration > 60 || isnan(SineSweep.Duration)
          error('Valor ingresado no válido.')
       end
    catch
        handles.SS.D.String = '30';
        SineSweep.Duration = str2double(get(handles.SS.D,'String'));
    end
    
 %% Errores de las frecuencias
         SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
         SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));
    try 
        if SineSweep.LowFreq<20||isnan(SineSweep.LowFreq)||SineSweep.LowFreq>=SineSweep.HighFreq||SineSweep.HighFreq>20000||isnan(SineSweep.HighFreq)
          error('Se esperan valores de frecuencia entre 20 y 20000 Hz. La frecuencia inferior no debe ser mayor que la superior.');
        end
    catch
        if SineSweep.LowFreq<20||isnan(SineSweep.LowFreq)
             handles.SS.F1.String = '88';
             SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
            
        end
        if SineSweep.HighFreq>20000||isnan(SineSweep.HighFreq)
            handles.SS.F2.String = '11314';
            SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));

        end
        if SineSweep.LowFreq>=SineSweep.HighFreq
            handles.SS.F1.String = '88';
            SineSweep.LowFreq = str2double(get(handles.SS.F1,'String'));
            handles.SS.F2.String = '11314';
            SineSweep.HighFreq = str2double(get(handles.SS.F2,'String'));
        end
    end
        [SineSweep] = SineSweep_FiltroInverso(SineSweep);
        setappdata(0,'SineSweep',SineSweep);
        handles.SScalc.txt4.String = 'El sine sweep y el filtro inverso han sido generados.';
        handles.SScalc.txt4.Visible = 'on';
end

function guardartable(~,~,handles)
% Funcion para guardar la tabla de resultados
% Se abren dos ventanas para guardar los resultados obtenidos con filtro de
% octavas y filtro de tercio de octavas por separado.
% Permite elegir el nombre y el formato del archivo.
% Por defecto propone crear un archivo excel llamado resultados

    try
        v = getappdata(0,'v');
        a = isempty(v);
        if a == 1
            error('Nooooooooooooooooooooooooooo');
        end
        filter = {'*.csv';'*.xlsx';'*.txt'};
        
        [file,path] = uiputfile(filter,'Guardar resultados con filtro de octavas','resultados_oct.xlsx');
        filename = fullfile(path,file);
        writetable(v.T1,filename,'WriteRowNames',true);
        
        [file,path] = uiputfile(filter,'Guardar resultados con filtro de tercio de octavas','resultados_terc.xlsx');
        filename = fullfile(path,file);
        writetable(v.T3,filename,'WriteRowNames',true);
    catch
        handles.ventanaerror.Visible = 'on';
        pause(1)
        set(handles.ventanaerror,'Visible','off');
    end
end

% Funciones para guardar el impulso generado
% Si se abre un archivo de una medición se guarda un impulso
% automaticamente con el nombre 'impulso_1'
% Si se abren multiples archivos de mediciones se guardan automaticamente todos los
% impulsos generados a partir de las mismas bajo el nombre de 'impulso_n'
% donde n representa el número del impulso guardado.

function guardarimpulso(~,~,handles)

    try
    matrix = getappdata(0,'matrix');
    a = isempty(matrix);
    if a == 1
        error('Nooooooooooooooooooooooooooo');
    end
    SaveIR(matrix);
    catch
        handles.ventanaerror.Visible = 'on';
        pause(1)
        set(handles.ventanaerror,'Visible','off');
    end
end

function SaveIR(matrix)

    for i = 1:matrix.n
        matrix.i = i;
        setappdata(0,'matrix',matrix);
        structfun(@write,matrix.RI,'UniformOutput',0);
    end
end

function write(h)

    matrix = getappdata(0,'matrix');
    filename = sprintf('Impulso_%d.wav',matrix.i);
    audiowrite(filename,h,matrix.SamplingFreq);
    
end
