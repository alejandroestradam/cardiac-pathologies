%%Acondicionamiento de biose�ales y bioim�genes
%%Proyecto final                     06/12/2019

close all; clc; 

%% Declaraci�n de variables
dentro_programa=1;
fs = 360;
T = 10;
t = 1/fs:1/fs:T;
ecg1=importdata('ECG1.mat'); %importar datos Ecg normal
ecg2=importdata('ECG2.mat'); %ECG complejo ventricular prematuro
ecg3=importdata('ECG3.mat'); %ECG taquicardia
ecg4=importdata('ECG4.mat'); %ECG fibrilaci�n atrial
ecga=ecg1(1:T*fs);
ecgb=ecg2(1:T*fs);
ecgc=ecg3(1:T*fs);
ecgd=ecg4(1:T*fs);

%% Filtrar registros de base de datos
dwtmode('per'); %Periodicidad: potencia 2
clc;
[b,a]=butter(2,[1.5 15]/(fs/2));
[d,c]=butter(2,150/(fs/2),'low');
ecg1_bs_butter=filter(b,a,ecga);
ecg1_bs_wavelet=wden(ecg1_bs_butter,'sqtwolog','s','sln',6,'db8');

ecg2_bs_butter=filter(d,c,ecgb);
ecg2_bs_wavelet=wden(ecg2_bs_butter,'sqtwolog','s','sln',8,'db8');

ecg3_bs_butter=filter(b,a,ecgc);
ecg3_bs_wavelet=wden(ecg3_bs_butter,'sqtwolog','s','sln',6,'db8');

ecg4_bs_butter=filter(b,a,ecgd);
ecg4_bs_wavelet=wden(ecg4_bs_butter,'sqtwolog','s','sln',6,'db8');

%% Graficar se�ales disponibles en base de datos sin filtrar y filtradas
figure(1)
subplot(2,2,1)
plot(t,ecga-950,'k')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG sano - Sin filtrar')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,2)
plot(t,ecgb-950,'r')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG CVP - Sin filtrar')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,3)
plot(t,ecgc-950,'b')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG Taquicardia - Sin filtrar')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,4)
plot(t,ecgd-950,'g')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG FA - Sin filtrar')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;

figure(2)
subplot(2,2,1)
plot(t,ecg1_bs_wavelet-950,'k')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG sano - Filtrado')
legend('bd8, nivel 6, sqtwolog soft')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,2)
plot(t,ecg2_bs_wavelet-950,'r')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG CVP - Filtrado')
legend('bd8, nivel 8, sqtwolog soft')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,3)
plot(t,ecg3_bs_wavelet-950,'b')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG Taquicardia - Filtrado')
legend('bd8, nivel 6, sqtwolog soft')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;
subplot(2,2,4)
plot(t,ecg4_bs_wavelet-950,'g')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG FA - Filtrado')
legend('bd8, nivel 6, sqtwolog soft')
xlim([1 10])
set(gca,'yticklabel',' ')
grid on;

%% Agregar ruido
ruido=wgn(1,3600,23);
ruido1=sin(2*pi*t);
added1=ecga+ruido1;
noised1=added1+ruido;
added1=ecgb+ruido1;
noised2=added1+ruido;
added1=ecgc+ruido1;
noised3=added1+ruido;
added1=ecgd+ruido1;
noised4=added1+ruido;

% Registros con ruido
figure(3)
subplot(2,2,1)
plot(t,noised1-950,'k')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG Sano - con ruido')
ylim([-50 450])
xlim([1 10])
set(gca,'yticklabel',' ')
grid on
subplot(2,2,2)
plot(t,noised2-950,'r')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG CVP - con ruido')
ylim([-350 550])
xlim([1 10])
set(gca,'yticklabel',' ')
grid on
subplot(2,2,3)
plot(t,noised3-950,'b')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG Taquicardia - con ruido')
ylim([-50 450])
xlim([1 10])
set(gca,'yticklabel',' ')
grid on
subplot(2,2,4)
plot(t,noised4-950,'g')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('ECG FA - con ruido')
ylim([-80 420])
xlim([1 10])
set(gca,'yticklabel',' ')
grid on

%%
while dentro_programa==1
    disp('BASE DE DATOS PARA DIAGN�STICO DE PATOLOG�AS DE ECG')
    disp('Seleccione archivo de registro de ECG: ')
    disp(' ')
    [name,dir]=uigetfile('.mat');
    archivo=strcat(dir,name);
    registro=importdata(archivo);
    ecg_usuario=registro(1:T*fs);

    %% Correlaci�n de ECG
    ecg_bd=[ecga; ecgb; ecgc; ecgd];
    ecg_bd_f=[ecg1_bs_wavelet; ecg2_bs_wavelet; ecg3_bs_wavelet; ecg4_bs_wavelet];
    for i=1:4
        C_ecg=corrcoef(ecg_usuario,ecg_bd_f(i,:)); %Matriz de correlaci�n
        cc_ecg(:,i)=abs(C_ecg(1,2)); %#ok
    end
    [val,pos]=max(cc_ecg);
    if val<=.1
        pos=0;
    end

    %% Diagn�stico de patolog�a
    switch pos
        case 1
            disp('No se detect� ninguna patolog�a, ECG sano.')
            disp(' ')
            diagn='ECG sano';
        case 2
            disp('Patolog�a diagnosticada: Complejo Ventricular Prematuro.')
            disp(' ')
            diagn='Patolog�a: Complejo Ventricular Prematuro';
        case 3
            disp('Patolog�a diagnosticada: Taquicardia.')
            disp(' ')
            diagn='Patolog�a: Taquicardia';
        case 4
            disp('Patolog�a diagnosticada: Fibrilaci�n Atrial.')
            disp(' ')
            diagn='Patolog�a: Fibrilaci�n Atrial';
        otherwise
            disp('Patolog�a no identificada.')
            disp(' ')
            disp('Fin del programa')
            dentro_programa=0;
    end
    
    %% Graficar Patolog�a diagnosticada y registro de usuario
    figure(4)
    subplot(2,1,1)
    plot(t,ecg_bd_f(pos,:)-950,'k')
    xlabel('Tiempo (s)')
    ylabel('Amplitud')
    title(diagn)
    xlim([1 10])
    set(gca,'yticklabel',' ')
    grid on;
    subplot(2,1,2)
    plot(t,ecg_usuario-950,'r')
    xlabel('Tiempo (s)')
    ylabel('Amplitud')
    title('ECG - Registro de usuario')
    xlim([1 10])
    set(gca,'yticklabel',' ')
    grid on;

    %% Procesar se�al
    disp('�Desea acondicionar la se�al?')
    disp('1. Si')
    disp('2. No')
    opc=input('Introduzca opci�n deseada: ');
    switch opc
            case 1
                dwtmode('per'); %Periodicidad: potencia 2
                [d,c]=butter(2,150/(fs/2),'low');
                switch pos
                    case 1
                        ecg_butter=filter(d,c,ecg_usuario);
                        ecg_wavelet=wden(ecg_butter,'sqtwolog','s','sln',10,'db8');
                    case 2
                        ecg_butter=filter(d,c,ecg_usuario);
                        ecg_wavelet=wden(ecg_butter,'sqtwolog','s','sln',10,'db8');
                    case 3
                        ecg_butter=filter(d,c,ecg_usuario);
                        ecg_wavelet=wden(ecg_butter,'sqtwolog','s','sln',10,'db8');
                    case 4
                        ecg_butter=filter(d,c,ecg_usuario);
                        ecg_wavelet=wden(ecg_butter,'sqtwolog','s','sln',10,'db8');
                end
                %% Graficar Patolog�a diagnosticada y registro de usuario
                figure(5)
                subplot(2,1,1)
                plot(t,ecg_usuario-950,'r')
                xlabel('Tiempo (s)')
                ylabel('Amplitud')
                title('ECG - Registro de usuario sin filtrar')
                xlim([1 10])
                set(gca,'yticklabel',' ')
                grid on;
                subplot(2,1,2)
                plot(t,ecg_wavelet-950,'b')
                xlabel('Tiempo (s)')
                ylabel('Amplitud')
                title('ECG - Registro de usuario filtrado')
                xlim([1 10])
                set(gca,'yticklabel',' ')
                grid on;
                dentro_programa=0;
            case 2
                disp(' ')
                disp('Fin del programa')
                dentro_programa=0;
    end
end