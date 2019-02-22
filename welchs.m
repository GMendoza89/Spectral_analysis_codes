%         Tarea Analisis espectral
%    Autor Gustavo David Mendoza pinto
%    fecha febrero 15 de 2019
%    
% 
clear
Fs = 1000;                             % Frecuencia de Muestreo          
N = 512;                              % Número de Muestras
F =300;                                % Frecuencia de la señal
F1 =370;                               % Frecuencia de la señal 2
M = 4096;                              % Muestras de la FFT
P = 50;                                % Número de repeticiones
n = zeros(1,N);                        % Vector de tiempo discreto
n(1) = 0;

K = 32;
L = N/K;

K1 = 16;
L1 = N/K1;


K2 = 8;
L2 = N/K1;

NP = .5;

w = hamming(L);
w1 = hamming(L1);
w2 = hamming(L2);

noverlap = L*NP;
noverlap1 = L1*NP;
noverlap2 = L2*NP;



for i = 2:N                            % Pasos de tiempo discreto
    n(i) = (1/Fs)*(i-1);    
end

for i=1 : P                            % Repeticiones del experimento
   X(:,i) = sin(2*pi*F*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + 2*sin(2*pi*F1*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + randn(N,1);
end

for j=1:P
PXX(:,j) = pwelch(X(:,j),w,noverlap,M);
end
PXXM=mean(PXX');

for j=1:P
PXX1(:,j) = pwelch(X(:,j),w1,noverlap1,M);
end
PXX1M=mean(PXX1');

for j=1:P
PXX2(:,j) = pwelch(X(:,j),w2,noverlap2,M);
end
PXX2M=mean(PXX2');
figure(1)                              % Gráfica 1 realicación de señal en un intervalo de tiempo
stem(X(1:200,1));

figure(2)
% Gráficas del análisis espectral
subplot(321)
plot(linspace(0,Fs/2,M/2),10*log10(PXX(1:M/2,:)));
grid on
title('wlch - 50 Repeticiones, 512 datos ventana rectangular de 16');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(322)
plot(linspace(0,Fs/2,M/2),10*log10(PXXM(1:M/2)));
plot(PXX2);
grid on
title('Welch Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');

subplot(323)
plot(linspace(0,Fs/2,M/2),10*log10(PXX1(1:M/2,:)));
grid on
title('wlch - 50 Repeticiones, 512 datos ventana rectangular de 32');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(324)
plot(linspace(0,Fs/2,M/2),10*log10(PXX1M(1:M/2)));
plot(PXX2);
grid on
title('Welch Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');

subplot(325)
plot(linspace(0,Fs/2,M/2),10*log10(PXX2(1:M/2,:)));
grid on
title('wlch - 50 Repeticiones, 512 datos ventana rectangular de 64');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(326)
plot(linspace(0,Fs/2,M/2),10*log10(PXX2M(1:M/2)));
plot(PXX2);
grid on
title('Welch Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');
set(2,'Color','White');
