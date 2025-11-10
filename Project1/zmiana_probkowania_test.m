%% resample_sine_script.m
% Parametry wejściowe
f0      = 150;        % [Hz] częstotliwość sinusa
Fs_in   = 25000;      % [Hz] oryginalna cz. próbkowania
Fs_out  = 100000;      % [Hz] docelowa cz. próbkowania
dur     = 0.0075;        % [s] czas trwania
A       = 1.0;        % amplituda
phi     = 0;          % [rad] faza początkowa

% --- Generacja sygnału wejściowego ---
t_in = (0:1/Fs_in:dur).';
x_in = A * sin(2*pi*f0*t_in + phi);

% --- Resampling / Interpolacja do Fs_out ---
if exist('resample','file') == 2
    % Wersja z filtracją (lepsza przy downsamplingu)
    [p,q] = rat(Fs_out/Fs_in, 1e-12);   % przybliżenie ułamkiem wymiernym
    x_out = resample(x_in, p, q);       % filtracja + zmiana Fs
    t_out = (0:numel(x_out)-1).' / Fs_out;
else
    % Fallback: czysta interpolacja czasowa (bez antyaliasingu)
    t_out = (0:1/Fs_out:dur).';
    % "spline" = gładko, "linear" = szybciej
    x_out = interp1(t_in, x_in, t_out, 'linear');
end

% --- Ostrzeżenie o aliasingu (gdy Fs_out za małe) ---
if f0 >= 0.5*Fs_out
    warning('Uwaga: f0 = %.1f Hz >= Fs_out/2 = %.1f Hz. Nastąpi aliasing.', f0, Fs_out/2);
end

% --- Porównanie z "idealnym" sinusem w nowej siatce ---
x_ref = A * sin(2*pi*f0*t_out + phi);
rmse  = sqrt(mean((x_out - x_ref).^2));
fprintf('RMSE względem idealnego sinusa w nowej siatce: %.3e\n', rmse);

% --- Podgląd ---
figure;
plot(t_in, x_in, 'DisplayName','Wejściowy (Fs_{in})'); hold on; grid on;
%stem(t_out, x_out, 'Marker','.', 'DisplayName','Po zmianie Fs','LineStyle','none');
plot(t_out, x_out, 'DisplayName','Po zmianie Fs','LineStyle','none');
xlabel('Czas [s]'); ylabel('Amplituda'); legend('Location','best');
title(sprintf('Resampling: f_0 = %.1f Hz, Fs_{in} = %d Hz -> Fs_{out} = %d Hz', f0, Fs_in, Fs_out));

% (opcjonalnie) zapis do WAV:
% audiowrite('sin_resampled.wav', x_out, Fs_out);
