clc
clear
close all

csv_read = readmatrix(uigetfile({'*.csv','CSV (*.CSV, )'}, 'Wybierz plik z danymi'));

data_read = csv_read(:,2);
fp = 100e3;
dt = 1/fp;
time = (dt : dt : size(data_read)/fp)';

X = fft(data_read);
X = X/max(X);

N = length(data_read);

figure
plot((fp*fp)/N*time, abs(X));
title('FFT sygnału');
xlim([0 fp/2]);
xlabel('Czestotliwosc [Hz]');
grid on;


figure
plot(time, data_read);
title('Podglad sygnału');
grid on;
xlabel('Czas [s]');
ylabel('Amplituda');
ylim([-1.1 1.1]);

%sound(data_read, fp);