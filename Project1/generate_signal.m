clc
clear
close all

fs = 100e3; %czestotiwosc probkowania
f_signal = 1e4; %czestotliwosc sygnalu
duration = 50; %czas trwania sygnalu w ms
duration = duration/1000;

ile_probek = 50; %przez jaki czas ma byc sygnal, reszta zera 
powtorzenia = 10; %ile razy powtorzyc sygnal

A = 1; %amplituda
t = (0 : 1/fs : duration-(1/fs))'; %tablica czasu dla sygnalu
y = A * sin(2*pi*f_signal*t); %generowanie sygnalu

ile_probek = ile_probek*100;
data = table(t, y);
temp_time = t(1:ile_probek);
temp_signal = A*sin(2*pi*f_signal*temp_time);
zero_time = t((ile_probek+1):end);
zero_signal = 0*sin(2*pi*f_signal*zero_time);

new_time = t;
new_signal = [temp_signal; zero_signal];

repeat_time = (0 : 1/fs : powtorzenia*(duration)-(1/fs))'; %rozszerz czas 10 razy

repeat_signal = [new_signal; new_signal];
for i = 1:(powtorzenia-2) %powtorz sygnal 10 razy
    repeat_signal = [repeat_signal; new_signal];
end

figure
plot(repeat_time, repeat_signal);
title('Przebieg sygnalu');
grid on;
xlabel('Czas [s]');
ylabel('Amplituda');
ylim([-1.1 1.1]);

X = fft(repeat_signal);
X = X/max(X);
N = length(repeat_signal);

figure
plot((fs*fs)/N*repeat_time, abs(X));
title('FFT sygnału');
xlim([0 fs/2]);
xlabel('Czestotliwosc [Hz]');
grid on;

Zapis = table(repeat_time, repeat_signal);

writetable(Zapis, 'new_signal.csv');