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

K = 4;
L = N/K;

K1 = 8;
L1 = N/K1;

for i = 2:N                            % Pasos de tiempo discreto
    n(i) = (1/Fs)*(i-1);    
end

for i=1 : P                            % Repeticiones del experimento
   X(:,i) = sin(2*pi*F*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + 2*sin(2*pi*F1*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + randn(N,1);
end

figure(1)                              % Gráfica 1 realicación de señal en un intervalo de tiempo
stem(X(1:200,1));



S = abs(fft(X,M)).^2/N;                % Cálculo del periodograma
SM=mean(S');                           % Promediados de la media de los periodogramas
set(1,'Color','White');


S1 = zeros(M,P);
for j=1:P
    for i = 1:L:(N-L)
        S1(:,j) = (S1(:,j) + abs(fft(X(i:(i+L),j),M)).^2/L)/K;                % Cálculo del periodograma
    end
end
SM1=mean(S1');                    % Promediados de la media de los periodogramas
set(1,'Color','White');

S2 = zeros(M,P);
for j=1:P
    for i = 1:L1:(N-L1)
        S2(:,j) = (S2(:,j) + abs(fft(X(i:(i+L1),j),M)).^2/L)/K1;                % Cálculo del periodograma
    end
end
SM2=mean(S2');                    % Promediados de la media de los periodogramas
set(1,'Color','White');


figure(2)                              % Gráficas del análisis espectral
subplot(321)
plot(linspace(0,Fs/2,M/2),10*log10(S(1:M/2,:)));
grid on
title('50 periodogramas con 5000 muestras');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(322)
plot(linspace(0,Fs/2,M/2),10*log10(SM(1:M/2)));
grid on
title('Periodograma Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');
set(2,'Color','White');

                            % Gráficas del análisis espectral
subplot(323)
plot(linspace(0,Fs/2,M/2),10*log10(S1(1:M/2,:)));
grid on
title('50 An‡lisis de Barllet con 512 muestras y 4 bloques');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(324)
plot(linspace(0,Fs/2,M/2),10*log10(SM1(1:M/2)));
grid on
title('An‡lisis de Barllet Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');

subplot(325)
plot(linspace(0,Fs/2,M/2),10*log10(S1(1:M/2,:)));
grid on
title('50 An‡lisis de Barllet con 512 muestras y 8 bloques');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(326)
plot(linspace(0,Fs/2,M/2),10*log10(SM2(1:M/2)));
grid on
title('An‡lisis de Barllet Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');

