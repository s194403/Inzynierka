% analiza_csv_fft.m
% Prosty odczyt CSV (t, x), wykres przebiegu i jego jednosided FFT

% === USTAWIENIA ===
%% --- Wybór pliku
[fn, fp] = uigetfile({'*.csv;*.txt','CSV/TXT (*.csv, *.txt)'}, 'Wybierz plik z danymi [czas, próbka]');
if isequal(fn,0), error('Przerwano wybór pliku.'); end
file = fullfile(fp, fn);
%filename = 'dane.csv';   % <-- podaj swoją nazwę pliku

% === WCZYTANIE DANYCH ===
M = readmatrix(file);           % auto-wykrywa separator (comma/semicolon)
if isempty(M) || size(M,2) < 2
    error('Plik "%s" nie zawiera co najmniej 2 kolumn z danymi.', filename);
end
t = M(:,1);                         % kolumna czasu (sekundy)
x = M(:,2);                         % kolumna wartości próbki

% Usunięcie ewentualnych NaN
valid = ~isnan(t) & ~isnan(x);
t = t(valid);
x = x(valid);

if numel(t) < 2
    error('Za mało danych: wczytano mniej niż 2 próbki.');
end

% (opcjonalnie) usunięcie składowej stałej
%x = detrend(x, 0);  % odejmuje średnią

% === SPRAWDZENIE KROKU PRÓBKOWANIA ===
dt = diff(t);
dt_med = median(dt);
Fs = 1/dt_med;                       % przybliżona częstotliwość próbkowania
if std(dt)/dt_med > 1e-3
    warning('Próbkowanie niejednostajne (~%.2f%% odchyłki). FFT traktuje dane jako równomierne.', ...
            100*std(dt)/dt_med);
end

% === WYKRES W DZIEDZINIE CZASU ===
figure; 
plot(t, x, 'LineWidth', 1);
grid on;
xlabel('Czas [s]');
ylabel('Amplituda [arb.]');
title('Przebieg w czasie');

% === FFT (jednostronne widmo amplitudowe) ===
N = numel(x);

% (opcjonalnie) okno Hann, ogranicza wyciek widma (lekkie zniekształcenie amplitudy)
% w = hann(N);
% xw = x .* w;
% Użyj bez okna dla najprostszej wersji:
xw = x;

X = fft(xw);
P2 = abs(X / N);                     % dwustronne
% części jednostronnej:
if mod(N,2)==0
    % N parzyste
    P1 = P2(1:N/2+1);
else
    % N nieparzyste
    P1 = P2(1:(N+1)/2);
end
P1(2:end-1) = 2*P1(2:end-1);         % skala jednostronna

% Oś częstotliwości (0 .. Fs/2)
f = Fs * (0:numel(P1)-1) / N;

% Wykres widma
figure;
plot(f, P1, 'LineWidth', 1);
grid on;
xlim([0, 1000]);
xlabel('Częstotliwość [Hz]');
ylabel('|X(f)|');
title('Jednostronne widmo amplitudowe (FFT)');
