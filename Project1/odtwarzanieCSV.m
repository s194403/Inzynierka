% analyze_csv_fft.m
% Odczyt: CSV [czas[s], sygnał], wykres w czasie i FFT (jednostronne widmo)

clear; clc;

%% --- Wybór pliku
[fn, fp] = uigetfile({'*.csv;*.txt','CSV/TXT (*.csv, *.txt)'}, 'Wybierz plik z danymi [czas, próbka]');
if isequal(fn,0), error('Przerwano wybór pliku.'); end
file = fullfile(fp, fn);

%% --- Wczytanie danych
% readmatrix poradzi sobie z nagłówkami — nie-numeryczne komórki -> NaN
raw = readmatrix(file);

% Usunięcie wierszy niepełnych / z NaN
raw = raw(all(~isnan(raw),2), :);

if size(raw,2) < 2
    error('Plik musi mieć co najmniej 2 kolumny: czas i wartość próbki.');
end

t = raw(:,1);
x = raw(:,2);

% Uporządkowanie po czasie (na wypadek permutacji) i usunięcie duplikatów czasu
[t, idx] = sort(t(:));
x = x(idx);
[tu, ia] = unique(t, 'stable');
xu = x(ia); 
t  = tu; 
x  = xu;

if numel(t) < 4
    error('Za mało próbek do analizy (min. 4).');
end

%% --- Sprawdzenie równomierności próbkowania
dt      = diff(t);
dt_med  = median(dt);
jitter  = max(abs(dt - dt_med)) / dt_med;  % względne odchylenie kroku
tol     = 1e-3;  % 0.1% tolerancji na nieregularność

is_uniform = jitter < tol;

if is_uniform
    fs_est = 1/dt_med;
    t_eq   = t;
    x_eq   = x;
else
    % Resampling do równomiernej siatki o kroku = mediana(diff(t))
    fs_est = 1/dt_med;
    t_eq   = (t(1):dt_med:t(end)).';
    % Interpolacja liniowa (możesz zmienić na 'pchip' dla gładszego przebiegu)
    x_eq   = interp1(t, x, t_eq, 'linear', 'extrap');
end

N = numel(x_eq);
dur = t_eq(end) - t_eq(1);

fprintf('Plik: %s\n', fn);
fprintf('Próbek (po obróbce): %d | czas: %.6f s | fs~: %.6f Hz | jitter: %.3g (%s)\n', ...
    N, dur, fs_est, jitter, ternary(is_uniform,'równomierne','NIERÓWNOMIERNE + resampling'));

%% --- (opcjonalnie) detrend/okno
do_detrend = true;
if do_detrend
    xw = detrend(x_eq, 'linear');  % usuń składową wolnozmienną/DC
else
    xw = x_eq;
end

% Proste okno Hann (poprawia estymację widma)
w  = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
xw = xw .* w;

%% --- FFT (jednostronne widmo amplitudowe)
Nfft = 2^nextpow2(N);
X    = fft(xw, Nfft);

% Normalizacja amplitudy: uwzględnij sumę okna, aby skala była bliższa realnej
win_gain = sum(w)/N;                      % średnia wagi okna
P2 = abs(X)/(N * win_gain);               % widmo dwustronne (amplituda)
P1 = P2(1:Nfft/2+1);
P1(2:end-1) = 2*P1(2:end-1);             % jednostronne
fax = fs_est*(0:(Nfft/2))/Nfft;

%% --- Wykresy
figure('Name','Analiza CSV: czas i FFT','Color','w');
tiledlayout(2,1,'TileSpacing','compact');

% 1) Czas (słupkowo przy dużej liczbie próbek przerzedzamy dla czytelności)
nexttile;
maxBars = 15000;
step    = max(1, ceil(N/maxBars));
stem(t_eq(1:step:end), x_eq(1:step:end), 'Marker','none'); grid on;
xlabel('Czas [s]'); ylabel('Amplituda');
title('Przebieg w dziedzinie czasu');

% 2) Widmo (amplituda liniowa; jeśli wolisz dB — patrz blok poniżej)
nexttile;
plot(fax, P1); grid on;
xlim([0, fs_est/2]);
xlabel('Częstotliwość [Hz]');
ylabel('|X(f)|');
title('Jednostronne widmo amplitudowe (FFT)');

%% --- (opcjonalnie) widmo w dBFS
%{
figure('Name','Widmo w dB','Color','w');
plot(fax, 20*log10(P1 + eps)); grid on;
xlim([0, fs_est/2]);
ylim([-160, 20]); % dostosuj wg poziomu sygnału
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dBFS]');
title('Jednostronne widmo amplitudowe (FFT) [dB]');
%}

%% --- Helper
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end