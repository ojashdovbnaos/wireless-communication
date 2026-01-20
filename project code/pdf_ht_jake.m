%% Combined SOS Channel Analysis: Amplitude, Phase & Time Variation (Single Figure)
clear; clc; close all;

%% 1) Simulation Parameters
N        = 128;            % number of sinusoids
fd       = 70;             % maximum Doppler shift [Hz]
fs       = 5000;           % sampling rate [Hz]
T        = 10;             % total simulation time [s]
t        = (0:1/fs:T-1/fs)'; % time vector [M×1], where M = fs*T

%% 2) Generate Random SOS Channel Parameters
rng('shuffle');            % optional: seed for randomness
theta    = 2*pi*rand(N,1);               % arrival angles [0,2π)
phi      = -pi + 2*pi*rand(N,1);         % initial phases [-π,π)

% Doppler frequencies f_n [N×1]
f_n      = fd * cos(theta);

%% 3) Generate Complex Envelope h(t) via Vectorized Sum-of-Sinusoids
exponent = 2*pi*(t * f_n.') + phi.';  % size M×N
h_t      = sum(exp(1j * exponent), 2) / sqrt(N);  % M×1, normalized energy

%% 4) Amplitude & Phase Extraction
amplitude = abs(h_t);      % Rayleigh-like amplitude
phase     = angle(h_t);    % Uniform phase in [–π,π)

%% 5) Single Figure with 3 Subplots
figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);

% (a) Amplitude Distribution vs. Rayleigh PDF
subplot(3,1,1);
histogram(amplitude, 50, 'Normalization', 'pdf');
hold on;
sigma = 1/sqrt(2);
x_pdf = linspace(0, max(amplitude), 200);
rayleigh_pdf = (x_pdf./sigma^2).*exp(-x_pdf.^2/(2*sigma^2));
plot(x_pdf, rayleigh_pdf, 'r-', 'LineWidth',1.5);
hold off;
xlabel('Amplitude |h(t)|');
ylabel('PDF');
title('(a) Amplitude Distribution vs. Rayleigh PDF');
legend('Simulation','Rayleigh PDF','Location','northeast');
grid on;

% (b) Phase Distribution (Uniform)
subplot(3,1,2);
edges = linspace(-pi, pi, 51);  % 50 bins
histogram(phase, edges, 'Normalization', 'pdf');
xlabel('Phase ∠h(t) [rad]');
ylabel('PDF');
title('(b) Phase Distribution');
xlim([-pi, pi]);
ylim([0, 1/(2*pi)*1.5]);
yline(1/(2*pi), 'r--', 'LineWidth',1.5);
legend('Simulation','Uniform PDF','Location','northeast');
grid on;

% (c) Time-Varying Channel Envelope
subplot(3,1,3);
plot(t, amplitude, 'LineWidth',1);
xlabel('Time [s]');
ylabel('Amplitude |h(t)|');
title('(c) Time-Varying Channel Envelope');
xlim([0, T]);
grid on;
