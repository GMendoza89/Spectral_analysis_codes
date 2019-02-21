%         Tarea Analisis espectral
%    Autor Gustavo David Mendoza pinto
%    fecha febrero 15 de 2019
%    
% 
clear
Fs = 1000;                             % Frecuencia de Muestreo          
N = 512;                              % Número de Muestras
F =200;                                % Frecuencia de la señal
F1 =370;                               % Frecuencia de la señal 2
M = 4096;                              % Muestras de la FFT
P = 50;                                % Número de repeticiones
n = zeros(1,N);                        % Vector de tiempo discreto
n(1) = 0;

L = 21;
L1 = 41;
L2 = 81;

W = ones(L,1);
W1 = ones(L1,1);
W2 = ones(L2,1);
W3 = hamming(L2);


for i = 2:N                            % Pasos de tiempo discreto
    n(i) = (1/Fs)*(i-1);    
end

for i=1 : P                            % Repeticiones del experimento
   X(:,i) = sin(2*pi*F*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + 2*sin(2*pi*F1*n(1,:)' + 2*pi*(rand(1,1)-0.5)) + randn(N,1);
end



% S=abs(fft(R.*window,M)).^2;

S = zeros(M,P);
for j=1:P
    R = xcorr(X,j);
    for i = 1:L:(size(R)-L)
        S(:,j) = (S(:,j) + abs(fft(W.*R(i:(i+L-1),1),M)).^2);                % Cálculo de la densidad espectral de potencia
    end
end
SM=mean(S');                    % Promediados de la media de los periodogramas

S1 = zeros(M,P);
for j=1:P
    R = xcorr(X,j);
    for i = 1:L:(size(R)-L1)
        S1(:,j) = (S1(:,j) + abs(fft(W1.*R(i:(i+L1-1),1),M)).^2);                % Cálculo del periodograma
    end
end
SM1=mean(S1');                    % Promediados de la media de los periodogramas

S2 = zeros(M,P);
for j=1:P
    R = xcorr(X,j);
    for i = 1:L2:(size(R)-L2)
        S2(:,j) = (S1(:,j) + abs(fft(W2.*R(i:(i+L2-1),1),M)).^2);                % Cálculo del periodograma
    end
end
SM2=mean(S2');                    % Promediados de la media de los periodogramas


S3 = zeros(M,P);
for j=1:P
    R = xcorr(X,j);
    for i = 1:L2:(size(R)-L2)
        S3(:,j) = (S1(:,j) + abs(fft(W3.*R(i:(i+L2-1),1),M)).^2);                % Cálculo del periodograma
    end
end
SM3=mean(S3');                    % Promediados de la media de las dencidades expectrales de potencia

figure(1)                              % Gráfica 1 realicación de señal en un intervalo de tiempo
stem(X(1:200,1));
set(1,'Color','White');

figure(2)
subplot(421)
plot(linspace(0,Fs/2,M/2),10*log10(S(1:M/2,:)));
grid on
title('Blacman - Tukey 50 repeticiones con 512 muestras y y ventana de rectangular de 21 valores');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(422)
plot(linspace(0,Fs/2,M/2),10*log10(SM(1:M/2)));
grid on
title('Blacman - Tukey Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');


subplot(423)
plot(linspace(0,Fs/2,M/2),10*log10(S1(1:M/2,:)));
grid on
title('Blacman - Tukey 50 repeticiones con 512 muestras y y ventana de rectangular de 41 valores');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(424)
plot(linspace(0,Fs/2,M/2),10*log10(SM1(1:M/2)));
grid on
title('Blacman - Tukey Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');


subplot(425)
plot(linspace(0,Fs/2,M/2),10*log10(S2(1:M/2,:)));
grid on
title('Blacman - Tukey 50 repeticiones con 512 muestras y y ventana de rectangular de 81 valores');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(426)
plot(linspace(0,Fs/2,M/2),10*log10(SM2(1:M/2)));
grid on
title('Blacman - Tukey Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');


subplot(427)
plot(linspace(0,Fs/2,M/2),10*log10(S3(1:M/2,:)));
grid on
title('Blacman - Tukey 50 repeticiones con 512 muestras y y ventana de Haming de 81 valores');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(428)
plot(linspace(0,Fs/2,M/2),10*log10(SM3(1:M/2)));
grid on
title('Blacman - Tukey Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 5000 muestras');
set(2,'Color','White');
