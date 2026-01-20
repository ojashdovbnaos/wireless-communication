%% Combined MEDS SOS Channel Analysis: Amplitude, Phase & Time Envelope
clear; clc; close all;

%% 1) Simulation Parameters
N    = 128;            % number of sinusoids
fd   = 70;             % maximum Doppler shift [Hz]
fs   = 5000;           % sampling rate [Hz]
T    = 10;             % total simulation time [s]
t    = (0:1/fs:T-1/fs)'; % time vector [fs*T × 1]

%% 2) Generate MEDS Doppler Frequencies & Random Phases
alpha = acos(linspace(-1 + 1/N, 1 - 1/N, N))';  % avoid endpoints in acos
f_n   = fd * cos(alpha);                        % [N×1]
phi   = 2*pi * rand(N,1);                       % random initial phases

%% 3) Construct Complex Envelope h(t) (Vectorized)
exponent = 2*pi*(t * f_n.') + phi.';            % [fs*T × N]
h_t      = sum(exp(1j * exponent), 2) / sqrt(N);% [fs*T × 1]

%% 4) Extract Amplitude & Phase
amplitude = abs(h_t);      % Rayleigh-like envelope
phase     = angle(h_t);    % should be uniform in [–π, π)

%% 5) Single Figure with 3 Subplots
figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);

% (a) Amplitude Distribution vs. Rayleigh PDF
subplot(3,1,1);
histogram(amplitude, 50, 'Normalization','pdf');
hold on;
sigma   = 1/sqrt(2);
x_pdf   = linspace(0, max(amplitude), 200);
ray_pdf = (x_pdf./sigma^2) .* exp(-x_pdf.^2/(2*sigma^2));
plot(x_pdf, ray_pdf, 'r-', 'LineWidth',1.5);
hold off;
xlabel('Amplitude |h(t)|');
ylabel('PDF');
title('(a) MEDS Amplitude Distribution');
legend('Simulation','Rayleigh PDF','Location','northeast');
grid on;

% (b) Phase Distribution vs. Uniform PDF
subplot(3,1,2);
edges = linspace(-pi, pi, 51);
histogram(phase, edges, 'Normalization','pdf');
hold on;
y_unif = 1/(2*pi) * ones(size(edges));
plot(edges, y_unif, 'r--','LineWidth',1.5);
hold off;
xlabel('Phase ∠h(t) [rad]');
ylabel('PDF');
title('(b) MEDS Phase Distribution');
xlim([-pi, pi]);
ylim([0, 1/(2*pi)*1.5]);
legend('Simulation','Uniform PDF','Location','northeast');
grid on;

% (c) Time-Varying Channel Envelope
subplot(3,1,3);
plot(t, amplitude, 'LineWidth',1);
xlabel('Time [s]');
ylabel('Amplitude |h(t)|');
title('(c) MEDS Time-Varying Envelope');
xlim([0, T]);
grid on;
