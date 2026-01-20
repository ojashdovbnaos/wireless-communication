%% Combined Doppler PSD Comparison: Jakes, MED, MEDS vs. Theoretical
clear; clc; close all;

%% 1) Parameter Definitions
fd   = 70;              % maximum Doppler shift [Hz]
fs   = 5000;            % sampling rate [Hz]
T    = 10;              % total simulation time [s]
N    = 64;              % number of sinusoids
t    = (0:1/fs:T-1/fs)';% time vector [Nx1]

%% 2) Compute Doppler Frequencies & Phases for Each Method
% Jakes method (random arrival angles)
theta_j   = 2*pi*rand(N,1);
f_jakes   = fd * cos(theta_j);
phi_jakes = -pi + 2*pi*rand(N,1);

% MED (Method of Equal Distances)
theta_med   = 2*pi*(0:N-1)'/N;
f_med       = fd * cos(theta_med);
phi_med     = -pi + 2*pi*rand(N,1);

% MEDS (Exact Doppler Spread)
n_meds      = (1:N)';
% invert Jakes PSD CDF: F(f) = (1/pi)(arcsin(f/fd)+pi/2) ⇒ f_meds
f_meds      = fd * sin(pi*((n_meds-0.5)-N/2)/N);
phi_meds    = -pi + 2*pi*rand(N,1);

methods   = {'Jakes','MED','MEDS'};
f_list    = {f_jakes, f_med, f_meds};
phi_list  = {phi_jakes, phi_med, phi_meds};

%% 3) Welch PSD Settings
win    = hamming(256);
nover  = 128;
nfft   = 1024;

%% 4) Generate & Plot PSD for Each SOS Method
figure; hold on;
for k = 1:numel(methods)
    f_n   = f_list{k};    % Doppler freqs (Nx1)
    phi_n = phi_list{k};  % phases (Nx1)
    
    % Sum N sinusoids to get fading process h(t)
    % vectorized: each column is one sinusoid
    exponent = 2*pi*(t * f_n.') + phi_n.';
    h_t      = sum(exp(1j * exponent), 2) / sqrt(N);
    
    % Estimate PSD via Welch
    [Pxx, F] = pwelch(h_t, win, nover, nfft, fs, 'centered');
    plot(F, 10*log10(Pxx), 'LineWidth', 1.2);
end

%% 5) Plot Theoretical Jakes Spectrum
f_th = linspace(-fd, fd, 2000);
S_th = 1 ./ (pi*fd * sqrt(1 - (f_th/fd).^2));  % valid for |f_th| < fd
plot(f_th, 10*log10(S_th), 'k--', 'LineWidth', 1.5);

%% 6) Figure Formatting
xlim([-fd, fd]);
ylim([-30, 0]);
xlabel('Doppler Shift f (Hz)');
ylabel('PSD (dB/Hz)');
title('Doppler Power Spectral Density: SOS Methods vs. Theoretical');
legend([methods, {'Theoretical Jakes'}], 'Location', 'SouthWest');
grid on;
