%% Combined Simulation & Theory: SOS-QPSK, OFDM-QPSK and QPSK Theoretical BER
clear; clc; close all;

%% 1) Simulation Parameters
snr_dB     = 0:5:20;           % SNR grid in dB
snr_linear = 10.^(snr_dB/10);  % linear SNR

%% 2) SOS-QPSK over COST207 Fading Channel
% – channel & symbol parameters
fs       = 5e3;      % sampling rate [Hz]
Tsim     = 10;       % total time [s]
fd_max   = 70;       % max Doppler [Hz]
N_sos     = 50;      % number of sinusoids
Rs       = 200;      % symbol rate [symbols/s]
sps      = fs/Rs;    % samples per symbol
t_vec    = linspace(0, Tsim, fs*Tsim);

path_delays = [0, 0.2, 0.4, 0.6]*1e-6;        % [s]
path_powers = [1, 0.63, 0.1, 0.01];
path_powers = path_powers/sum(path_powers);

% – generate per‐symbol fading h_sym (1×Nsym)
h_total = zeros(numel(path_delays), numel(t_vec));
for p = 1:numel(path_delays)
    theta = 2*pi*rand(1,N_sos);
    phi   = -pi + 2*pi*rand(1,N_sos);
    h_p   = sum(exp(1j*(2*pi*fd_max*cos(theta).' * t_vec + phi.')), 1);
    h_total(p,:) = sqrt(path_powers(p)) * (h_p / sqrt(N_sos));
end
h_sym = sum(h_total(:,1:sps:end),1);
Nsym  = numel(h_sym);

% – QPSK symbol generation
bits_q  = randi([0 1], 1, 2*Nsym);
sym_idx = bi2de(reshape(bits_q,2,[])','left-msb');
sym_q   = pskmod(sym_idx, 4, pi/4, 'gray') .';  % 1×Nsym

% – simulate BER vs. SNR
ber_sos = zeros(size(snr_dB));
for k = 1:numel(snr_dB)
    rx        = h_sym .* sym_q;
    rx_noisy  = awgn(rx, snr_dB(k), 'measured');
    rx_equal  = rx_noisy ./ h_sym;
    idx_hat   = pskdemod(rx_equal, 4, pi/4, 'gray');
    bits_hat  = de2bi(idx_hat, 2, 'left-msb')';
    bits_hat  = bits_hat(:).';
    [~, ber_sos(k)] = biterr(bits_q, bits_hat);
end

%% 3) OFDM-QPSK over COST207 Rayleigh Channel
% – OFDM / channel parameters
Nfft        = 64;              % number of subcarriers
CP_len      = 16;              % cyclic prefix length
M           = 4;               % QPSK
bps         = log2(M);         % bits per symbol
Nsym_ofdm   = 1e3;             % OFDM symbols per SNR
pd          = [0,3e-6,5e-6,8e-6];
pg_dB       = [0,-3,-6,-9];
fd_chan     = 100;             % max Doppler [Hz]
Fs_chan     = 1e6;             % channel sampling rate [Hz]

raylChan = comm.RayleighChannel( ...
    'SampleRate',         Fs_chan, ...
    'PathDelays',         pd, ...
    'AveragePathGains',   pg_dB, ...
    'MaximumDopplerShift',fd_chan, ...
    'NormalizePathGains', true, ...
    'PathGainsOutputPort',true);

ber_ofdm = zeros(size(snr_dB));
totalBits = Nsym_ofdm * Nfft * bps;

for k = 1:numel(snr_dB)
    reset(raylChan);
    bit_errors = 0;
    
    for sym_i = 1:Nsym_ofdm
        % generate bits & QPSK
        bits_tx = randi([0 1], Nfft*bps, 1);
        idx_tx   = bi2de(reshape(bits_tx,bps,[])','left-msb');
        qpsk_tx  = pskmod(idx_tx, M, pi/4);
        
        % OFDM modulate
        tx_ifft = ifft(qpsk_tx, Nfft);
        tx_cp   = [tx_ifft(end-CP_len+1:end); tx_ifft];
        
        % pass through channel + AWGN
        [rx_cp, pathGains] = raylChan(tx_cp);
        rx_noisy          = awgn(rx_cp, snr_dB(k), 'measured');
        
        % OFDM demodulate
        rx_fft = fft(rx_noisy(CP_len+1:end), Nfft);
        
        % ideal ZF equalization
        gains = mean(pathGains(CP_len+1:end,:), 1);
        h_tap = zeros(Nfft,1);
        samp_delays = round(pd*Fs_chan)+1;
        for p = 1:numel(samp_delays)
            if samp_delays(p)<=Nfft
                h_tap(samp_delays(p)) = gains(p);
            end
        end
        H = fft(h_tap, Nfft);
        eq_sym = rx_fft ./ H;
        
        % demod & count errors
        idx_hat = pskdemod(eq_sym, M, pi/4);
        bits_hat = reshape(de2bi(idx_hat, bps, 'left-msb').',[],1);
        bit_errors = bit_errors + sum(bits_hat~=bits_tx);
    end
    
    ber_ofdm(k) = bit_errors / totalBits;
end

%% 4) Theoretical Curves
% convert to Eb/N0 assuming snr_dB = Es/N0
ebn0 = snr_linear/2;

% A) Rayleigh‐fading average QPSK BER
ber_theory_rayleigh = 0.5*(1 - sqrt( ebn0 ./ (1 + ebn0) ));

% B) AWGN QPSK BER  (for comparison)
ber_theory_awgn = qfunc( sqrt(2*ebn0) );

%% 5) Plot All Curves
figure;
semilogy(snr_dB, ber_sos,  'o-', 'LineWidth',1.4); hold on;
semilogy(snr_dB, ber_ofdm,'s--','LineWidth',1.4);

semilogy(snr_dB, ber_theory_awgn,     '^-','LineWidth',1.2);
grid on;

xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('SOS-QPSK vs OFDM-QPSK Simulation and QPSK Theory');
legend( ...
  'SOS-QPSK + ZF (sim)', ...
  'OFDM-QPSK + ZF (sim)', ...
  'AWGN QPSK Theory', ...
  'Location','southwest' );

axis([snr_dB(1), snr_dB(end), 1e-4, 1]);
