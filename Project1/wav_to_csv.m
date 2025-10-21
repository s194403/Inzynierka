clc
close all
clear

[data_read, fp] = audioread(uigetfile({'*.wav','WAV (*.wav, )'}, 'Wybierz plik z danymi'));

dt = 1/fp;

time = (dt : dt : size(data_read)/fp);
time_write = time';

data_write = (data_read(:,1) + data_read(:,2)) /2;

Table = table(time_write, data_write);

writetable(Table, 'sample.csv');

X = fft(data_write);

N = length(data_write);

figure
title('FFT sygnału');
plot((fp*fp)/N*time, abs(X));
xlim([0 2000]);
xlabel('Czestotliwosc [Hz]');
grid on;

