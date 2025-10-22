clc
clear
close all

csv_read = readmatrix(uigetfile({'*.csv','CSV (*.CSV, )'}, 'Wybierz plik z danymi'));

data_read = csv_read(:,2);
fp = 48000;
dt = 1/fp;
time = (dt : dt : size(data_read)/fp)';

X = fft(data_read);

N = length(data_read);

figure
title('FFT sygnału');
plot((fp*fp)/N*time, abs(X));
xlim([0 2000]);
xlabel('Czestotliwosc [Hz]');
grid on;


figure
title('Podglad sygnału');
plot(time, data_read);
grid on;
xlabel('Czas [s]');
ylabel('Amplituda');

sound(data_read, fp);