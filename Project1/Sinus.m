% gen_sine_csv.m
% Prosty generator sinusa + uśrednianie w oknach + zapis do CSV (bez WAV)
clear; clc;

%% Parametry
fs       = 10000;            % [Hz] cz. próbkowania
f        = 100;              % [Hz] cz. sinusa
duration = 50*(1/f);         % [s] czas trwania (np. 50 okresów)
A        = 0.8;              % amplituda (0..1)
phi      = 0;                % [rad] faza początkowa
stereo   = false;            % true = stereo (L=R), false = mono
basename = 'sine_100Hz';     % baza nazwy plików wynikowych (bez rozszerzeń)

%% Generacja
t = (0:1/fs:duration).';           % kolumna czasu
y = A * sin(2*pi*f*t + phi);       % sygnał mono

%% --- UŚREDNIANIE w stałych oknach ---
win_ms        = 1;                                  % długość okna [ms]
win_samples   = max(1, round((win_ms/1000)*fs));    % próbki/okno
N             = size(y,1);
C             = size(y,2);                          % 1 = mono, 2 = stereo
num_blocks    = ceil(N / win_samples);

y_avg = zeros(num_blocks, C);
starts = 1:win_samples:N;
ends   = min(starts + win_samples - 1, N);
for k = 1:numel(starts)
    y_avg(k, :) = mean(y(starts(k):ends(k), :), 1);
end

% Czas próbek uśrednionych = środki okien
t_avg = ((starts(:) + ends(:) - win_samples -1)/2) / fs;         % [s], kolumna

% „fs” uśrednionego strumienia (tylko informacyjnie)
fs_avg = fs / win_samples;     % bez zaokrąglania; realnie próbki co win_ms

%% --- ZAPIS do CSV (czas + próbki)
% 1) oryginał
T_orig = table(t, y(:,1), 'VariableNames', {'time_s','sample'});
writetable(T_orig, sprintf('%s_orig.csv', basename));


fprintf('Zapisano CSV: %s_orig.csv  (N=%d)\n', basename, height(T_orig));

%% --- FFT: porównanie widma oryginału i wersji uśrednionej (kanał 1)
ch = 1;
y0 = y(:, ch);
y1 = y_avg(:, ch);

% Oryginał
L0 = length(y0);
Y0 = fft(y0);
P20 = abs(Y0/L0);
P10 = P20(1:floor(L0/2)+1);
P10(2:end-1) = 2*P10(2:end-1);
f0 = fs*(0:floor(L0/2))/L0;

% Po uśrednianiu (uwaga: nierówne „fs” jeśli ostatnie okno krótsze;
% tu przyjmujemy nominalny fs_avg = fs/win_samples)
L1 = length(y1);
Y1 = fft(y1);
P21 = abs(Y1/L1);
P11 = P21(1:floor(L1/2)+1);
P11(2:end-1) = 2*P11(2:end-1);
f1 = fs_avg*(0:floor(L1/2))/L1;

%% --- Wykresy: czas i widmo
figure('Name','Czas','Color','w');
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(t, y(:,1)); grid on; xlabel('Czas [s]'); ylabel('Amplituda');
title('Oryginał');
nexttile;
plot(t_avg, y_avg(:,1)); grid on; xlabel('Czas [s]'); ylabel('Amplituda');
title(sprintf('Uśredniony (okno = %d ms)', win_ms));
%linkaxes(findall(gcf,'Type','axes'),'x');

figure('Name','Widmo','Color','w');
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(f0, P10); grid on; xlim([0, 2*f]); xlabel('Częstotliwość [Hz]'); ylabel('|Y(f)|');
title(sprintf('FFT oryginału (fs = %g Hz, f = %g Hz)', fs, f));
nexttile;
plot(f1, P11); grid on; xlim([0, 2*f]); xlabel('Częstotliwość [Hz]'); ylabel('|Y(f)|');
title(sprintf('FFT uśrednionego (fs_{avg} ≈ %g Hz, okno = %d ms)', fs_avg, win_ms));

% (opcjonalnie) tłumienność boxcar
T = win_ms/1000;
theo_att = abs(sin(pi*f*T)/(pi*f*T));   % |sinc(pi f T)|
fprintf('Teoretyczna tłumienność przy %g Hz i oknie %d ms: %.4f (%.1f%% amplitudy)\n', ...
        f, win_ms, theo_att, 100*theo_att);