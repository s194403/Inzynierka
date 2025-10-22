clc
close all
clear

file = uigetfile({'*.wav','WAV (*.wav, )'}, 'Wybierz plik z danymi');

[data_read, fp] = audioread(file);
info = audioinfo(file)

dt = 1/fp;

time = (dt : dt : size(data_read)/fp);
time_write = time';

data_write = (data_read(:,1)); % + data_read(:,2)) /2; %w zaleznosci czy mono czy stereo

Table = table(time_write, data_write);

writetable(Table, 'kryzys.csv');

X = fft(data_write);

N = length(data_write);

figure
title('FFT sygnału');
plot((fp*fp)/N*time, abs(X));
xlim([0 2000]);
xlabel('Czestotliwosc [Hz]');
grid on;
