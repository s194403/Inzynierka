clc
clear
close all

csv_read = readmatrix(uigetfile({'*.csv','CSV (*.CSV, )'}, 'Wybierz plik z danymi'));

data_read = csv_read(:,2);
fp = 44100;
dt = 1/fp;
time = (dt : dt : size(data_read)/fp)';

figure
title('Podglad sygnału');
plot(time, data_read);
grid on;
xlabel('Czas [s]');

sound(data_read, fp);