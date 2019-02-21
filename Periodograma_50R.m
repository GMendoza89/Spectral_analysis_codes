%% Periodograma con cincuenta repeticiones
clear all
N = 8192;
M = 16384;
P=50;
x = 5*randn(N,50);
figure(1)
S=abs(fft(x,M)).^2/N;
SM=mean(S');
subplot(211)
plot(linspace(0,1,M/2),10*log10(S(1:M/2,:)));
grid on
title('50 periodogramas con 8192 muestras');
ylabel('Densidad espectral de Potencia (dB)');
xlabel('Frecuencia');
subplot(212)
plot(linspace(0,1,M/2),10*log10(SM(1:M/2)))
grid on
title('Periodograma Promediado');
ylabel('Densidad Espectral de Potencia (dB)');
xlabel('50 realizaciones del periodograma de WGN de 8192 muestras');
set(1,'Color','White');
