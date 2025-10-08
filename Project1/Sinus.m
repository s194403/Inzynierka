% gen_sine_wav.m
% Prosty generator pliku WAV z sinusoidą + uśrednianie co 5 ms + FFT

%% Parametry
fs       = 48000;          % [Hz] częstotliwość próbkowania
f        = 50;             % [Hz] częstotliwość sinusa
duration = 0.5;            % [s] czas trwania
A        = 0.8;            % amplituda (0..1)
phi      = 0;              % [rad] faza początkowa
bits     = 16;             % rozdzielczość zapisu: 16/24/32
stereo   = false;          % true = stereo (L=R), false = mono
filename = 'sine_440Hz.wav'; % (nazwa może nie zgadzać się z f; tylko etykieta)

%% Generacja
t = (0:1/fs:duration).';           % kolumna czasu
y = A * sin(2*pi*f*t + phi);       % sygnał mono

% (opcjonalnie) stereo: sklonuj na 2 kanały
if stereo
    y = [y y];                      % [N x 2]
end

% Zabezpieczenie przed clippingiem
maxabs = max(abs(y),[],'all');
if maxabs > 1
    y = y / maxabs;
    warning('Sygnał znormalizowano, aby uniknąć przesterowania.');
end

%% --- UŚREDNIANIE CO 5 ms ---
win_ms        = 5;                                  % długość okna w milisekundach
win_samples   = max(1, round((win_ms/1000)*fs));    % liczba próbek w 5 ms
N             = size(y,1);
C             = size(y,2);                          % 1 = mono, 2 = stereo
num_blocks    = ceil(N / win_samples);

y_avg = zeros(num_blocks, C);
for k = 1:num_blocks
    s = (k-1)*win_samples + 1;
    e = min(k*win_samples, N);
    y_avg(k, :) = mean(y(s:e, :), 1);              % średnia w oknie dla każdego kanału
end

% Efektywna fs po uśrednianiu (1 próbka na okno 5 ms)
fs_avg = round(fs / win_samples);                   % ≈ 200 Hz przy fs=48 kHz

%% Zapis do WAV
% 1) oryginał
audiowrite(filename, y, fs, 'BitsPerSample', bits);

% 2) wersja uśredniona co 5 ms
[pth, nam, ~] = fileparts(filename);
filename_avg = fullfile(pth, [nam, '_avg_', num2str(win_ms), 'ms.wav']);
audiowrite(filename_avg, y_avg, fs_avg, 'BitsPerSample', bits);

%% --- FFT: porównanie widma oryginału i wersji uśrednionej ---
ch = 1;                         % który kanał analizować (1 = lewy/mono)
y0 = y(:, ch);
y1 = y_avg(:, ch);

% Oryginał
L0 = length(y0);
Y0 = fft(y0);
P20 = abs(Y0/L0);
P10 = P20(1:floor(L0/2)+1);
P10(2:end-1) = 2*P10(2:end-1);
f0 = fs*(0:floor(L0/2))/L0;

% Po uśrednianiu
L1 = length(y1);
Y1 = fft(y1);
P21 = abs(Y1/L1);
P11 = P21(1:floor(L1/2)+1);
P11(2:end-1) = 2*P11(2:end-1);
f1 = fs_avg*(0:floor(L1/2))/L1;

% Teoretyczna tłumienność (filtr prostokątny 5 ms)
T = win_ms/1000;
theo_att = abs(sin(pi*f*T)/(pi*f*T));   % |sinc(pi f T)|
fprintf('Teoretyczna tłumienność przy %d Hz i oknie %d ms: %.4f (%.1f%% amplitudy)\n', ...
        f, win_ms, theo_att, 100*theo_att);

% Wykresy FFT (skala liniowa)
figure('Name','Widmo amplitudowe');
tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(f0, P10); grid on;
xlim([0, min(fs/2, 100)]); % pokaż do 500 Hz (dla przejrzystości)
xlabel('Częstotliwość [Hz]');
ylabel('|Y(f)|');
title(sprintf('Oryginał (fs = %d Hz, f = %.2f Hz)', fs, f));

nexttile;
plot(f1, P11); grid on;
xlim([0, min(fs_avg/2, 100)]);
xlabel('Częstotliwość [Hz]');
ylabel('|Y(f)|');
title(sprintf('Po uśrednianiu (fs = %d Hz, okno = %d ms)', fs_avg, win_ms));

%% (opcjonalnie) odsłuch
% sound(y, fs);
% sound(y_avg, fs_avg);

fprintf('Zapisano oryginał:   %s  | fs=%d Hz, f=%.2f Hz, czas=%.2f s, %s, %d-bit\n', ...
    filename, fs, f, duration, ternary(stereo,'stereo','mono'), bits);
fprintf('Zapisano uśredniony: %s  | fs_avg=%d Hz (co %d ms), %s, %d-bit\n', ...
    filename_avg, fs_avg, win_ms, ternary(stereo,'stereo','mono'), bits);

% Po policzeniu y_avg:
win_s   = win_ms/1000;             % 0.005 s
starts  = 1:win_samples:N;         % indeksy początków okien
ends    = min(starts + win_samples - 1, N);   % indeksy końców okien
t_avg   = ((starts + ends)/2 - 1).' / fs;     % czasy środków okien [s] (kolumna)


figure;
subplot(2,1,1);
plot(t, y(:,1)); grid on;           % jeśli stereo, bierz kanał 1
title('Oryginał'); xlabel('Czas [s]'); ylabel('Amplituda');
subplot(2,1,2);
plot(t_avg, y_avg(:,1)); grid on;
title('Po uśrednianiu'); xlabel('Czas [s]'); ylabel('Amplituda');

% --- helper: operator warunkowy jak w C (tylko do ładnego printa)
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
