% analyze_wav_fft.m
% Wczytuje plik WAV, rysuje przebieg w czasie (słupkowo) i widmo FFT.

clear; clc;

%% --- Wybór pliku
[fn, fp] = uigetfile({'*.wav','WAV files (*.wav)'}, 'Wybierz plik WAV');
if isequal(fn,0)
    error('Przerwano wybór pliku.');
end
file = fullfile(fp, fn);

%% --- Wczytanie audio i metadanych
[x, fs] = audioread(file);        % x: [N x C]
info    = audioinfo(file);
N       = size(x,1);
C       = size(x,2);

fprintf('Plik: %s\nfs: %d Hz | kanały: %d | bit depth: %d | długość: %.2f s\n', ...
    fn, fs, C, info.BitsPerSample, N/fs);

%% --- Wybór kanału do analizy (1 dla mono/lewego)
chan = 1;
x_ch = x(:, min(chan, C));

%% --- Wykres w dziedzinie czasu (słupkowy)
t = (0:N-1).'/fs;

% Aby wykres był czytelny, ograniczamy liczbę „słupków”
maxBars = 8000;
step    = max(1, ceil(N/maxBars));
tb      = t(1:step:end);
xb      = x_ch(1:step:end);

%% --- FFT (jednostronne widmo amplitudowe)
% Okno Hann (bez zależności od toolboxów)
w  = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
xw = x_ch .* w;

Nfft = 2^nextpow2(N);             % zero-padding do potęgi 2 (ładniejsze widmo)
X    = fft(xw, Nfft);
P2   = abs(X/N);                   % widmo dwustronne (amplituda)
P1   = P2(1:Nfft/2+1);             % jednostronne
P1(2:end-1) = 2*P1(2:end-1);       % korekta energii
fax  = fs*(0:(Nfft/2))/Nfft;

%% --- Rysowanie
figure('Name','Analiza WAV','Color','w');
tiledlayout(2,1,'TileSpacing','compact');

% 1) Czas
nexttile;
%stem(tb, xb, 'Marker','none');
plot(tb, xb);
grid on;
xlabel('Czas [s]');
ylabel('Amplituda');
title(sprintf('Przebieg w czasie (kanał %d, %s) — %s', ...
    chan, ternary(C==1,'mono','stereo'), fn), 'Interpreter','none');

% 2) Widmo
nexttile;
plot(fax, P1); grid on;
xlim([0, 1000]);
xlabel('Częstotliwość [Hz]');
ylabel('|X(f)|');
title('Jednostronne widmo amplitudowe (FFT)');

%% --- (opcjonalnie) wersja dB zamiast liniowej amplitudy
%{
nexttile;
plot(fax, 20*log10(P1 + eps)); grid on;
xlim([0, fs/2]);
ylim([-160, 20]); % dostosuj do poziomów sygnału
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dBFS]');
title('Jednostronne widmo amplitudowe (FFT) [dB]');
%}

%% --- Helper (operator warunkowy jak w C)
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
