clc
close all
clear

[data_read, fp] = audioread(uigetfile({'*.wav','WAV (*.wav, )'}, 'Wybierz plik z danymi'));

time = (1/fp : 1/fp : size(data_read)/fp);
time = time';

data_write = (data_read(:,1) + data_read(:,2)) /2;

Table = table(time, data_write);

writetable(Table, 'sample.csv');

