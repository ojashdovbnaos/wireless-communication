%% compare_sos_tdl.m
% Compare TDL h(τ,t) responses for three SOS parameterization methods
clear; clc; close all;

%% 1) Simulation Parameters
N       = 128;                      % number of sinusoids
fd      = 70;                       % maximum Doppler shift [Hz]
fs      = 5000;                     % sampling rate [Hz]
T       = 10;                       % simulation duration [s]
t       = linspace(0, T, fs*T);     % time vector [1×L]
L       = numel(t);

% COST207 4-path delays & powers
taus    = [0, 0.2e-6, 0.4e-6, 0.6e-6]; % path delays [s]
powers  = [1.0, 0.63, 0.10, 0.01];
powers  = powers./sum(powers);         % normalized
P       = numel(taus);
taus_us = taus*1e6;                    % for y-axis [μs]

% decimation for plotting clarity
M   = 50;
idx = 1:M:L;
t2  = t(idx);

%% 2) Preallocate Channel Matrices
H_rand = zeros(P, L);
H_med  = zeros(P, L);
H_meds = zeros(P, L);

%% 3) Generate h(τ,t) for Each Method & Path
for p = 1:P
    w = sqrt(powers(p));  % path weight

    % 3.1) Random theta method
    theta = 2*pi*rand(N,1);
    phi   = -pi + 2*pi*rand(N,1);
    f_rand = fd*cos(theta);
    E_rand = 2*pi*(t.'*f_rand.') + phi.';       % [L×N]
    h_rand = sum(exp(1j*E_rand), 2).' / sqrt(N);% [1×L]
    H_rand(p,:) = w * h_rand;

    % 3.2) MED method
    theta_m = 2*pi*((1:N).'-0.5)/N;
    f_med   = fd*cos(theta_m);
    phi_m   = -pi + 2*pi*rand(N,1);
    E_med   = 2*pi*(t.'*f_med.') + phi_m.';
    h_med   = sum(exp(1j*E_med), 2).' / sqrt(N);
    H_med(p,:)  = w * h_med;

    % 3.3) MEDS method
    alpha  = acos(linspace(-1+1/N,1-1/N,N)).';
    f_meds = fd*cos(alpha);
    phi_s  = -pi + 2*pi*rand(N,1);
    E_meds = 2*pi*(t.'*f_meds.') + phi_s.';
    h_meds = sum(exp(1j*E_meds), 2).' / sqrt(N);
    H_meds(p,:) = w * h_meds;
end


%% 4) Plot with tiledlayout
methods = {'Random θ','MED','MEDS'};
H_list  = {H_rand, H_med, H_meds};

% 改成 1 行 3 列
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    ax = nexttile;
    hold(ax,'on');
    for p = 1:P
        plot3(ax, t2, taus_us(p)*ones(size(t2)), ...
              abs(H_list{k}(p,idx)), 'LineWidth',1.2);
    end
    view(ax,45,30);
    xlabel(ax,'Time (s)');
    ylabel(ax,'Delay τ (µs)');
    zlabel(ax,sprintf('|h_{%s}|', lower(methods{k})));
    title(ax,methods{k});
    grid(ax,'on');
end
