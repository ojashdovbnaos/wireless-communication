%% SOS‐TDL Channel Impulse Response Visualization
% Generate and plot h(τ,t) for a 4‐tap COST207 TDL using SOS

%% 1) Clear and Parameters
clear; close all; clc;

% SOS parameters
N       = 128;            % number of sinusoids per path
fd      = 70;             % max Doppler shift [Hz]
fs      = 5000;           % sampling rate [Hz]
T       = 10;             % duration [s]
t       = (0:1/fs:T-1/fs)';% time vector [L×1], L = fs*T

% COST207 4‐tap delays & normalized powers
taus    = [0, 0.2e-6, 0.4e-6, 0.6e-6]; % delays [s]
powers  = [1.0, 0.63, 0.10, 0.01];
powers  = powers / sum(powers);
P       = numel(taus);

% Decimation for plotting clarity
M       = 50;                    
idx     = 1:M:numel(t);         
t2      = t(idx);                

% Convert delays to µs for plotting
taus_us = taus * 1e6;            

colors  = lines(P);              % distinct colors

%% 2) Generate h(τ,t) Matrix
H = zeros(P, numel(t));          % preallocate [P × L]

for p = 1:P
    % random arrival angles and phases
    theta = 2*pi*rand(N,1);
    phi   = -pi + 2*pi*rand(N,1);
    
    % compute Doppler freqs and exponent matrix [L×N]
    f_n     = fd * cos(theta);
    exponent= 2*pi*(t * f_n.') + phi.';
    
    % SOS sum and power weighting
    h       = sum(exp(1j*exponent), 2) / sqrt(N);  % [L×1]
    H(p,:)  = sqrt(powers(p)) * h.';               % store row
end

%% 3) 3D Plot of |h(τ,t)|
figure('Name','h(\tau,t) — SOS+TDL Response','NumberTitle','off');
hold on;
for p = 1:P
    plot3(...
      t2, ...                         % time axis
      taus_us(p)*ones(size(t2)), ... % delay axis
      abs(H(p,idx)), ...             % magnitude
      'Color',   colors(p,:), ...
      'LineWidth',1.2 ...
    );
end
hold off;

% Formatting
xlabel('Time (s)');
ylabel('Delay \tau (µs)');
zlabel('|h(\tau,t)|');
title('Unit Impulse Response of COST207 Channel via SOS');
grid on;
view(45,30);

% Legend with delay labels
legend_labels = arrayfun(@(d) sprintf('%.1f µs', d), taus_us, 'UniformOutput', false);
legend(legend_labels, 'Location','northeast');
