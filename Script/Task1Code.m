clear all; close all; clc;

MC = 1000;   % number of Monte Carlo runs for data
MC_channel = 1000; % number of Monte Carlo runs for channel
N_TX = 2; % number of TX antennas
N_RX = 2; % number of RX antennas

EbN0dB = 0:0.1:10;  % SNR vector in dB
NoSNRPoints = length(EbN0dB);   % length of SNR vector
EbN0 = 10.^(EbN0dB./10);  % SNR in linear scale
C = ... ;          % Find noise variance 
sig = ... ;        % Find noise standard deviation 

A = sign(randn(N_TX,MC));  % random data for all runs, the same for all SNR points
C_A = ... ;                % Find data variance (variance of A)
ASE_LS = zeros(NoSNRPoints,MC_channel);  % initialize the LS average square error (ASE)
ASE_LMMSE = zeros(NoSNRPoints,MC_channel);  % LMMSEE ASE
MSE_LMMSE = zeros(NoSNRPoints,MC_channel);  % theoretical MSE for LMMSE

for i = 1:NoSNRPoints        % generate noise realization for SNR values using
    for c = 1:MC_channel     %generate channel realization 
        H = ...;              % generate complex channel response with total variance 1
        w = ...;              % generate complex noise with different sig for different SNR
        X = H*A+ w;
        LS = ...;             % Implement LS estimator and take only real part which contain the data
        LMMSE = ...;          % Implement LMMSE estimator and take only real part which contain the data
        ASE_LS(i,c) = mean((abs(LS-A)).^2,'all');  % Calculate average square error for LS
        ASE_LMMSE(i,c) = ...;                       % Calculate average square error for LMMSE
        MSE(i,c) = sum(abs(diag(real(...)))); % Implement theoritical MSE for LMMSE by taking real parts of diagonal elements
    end
end

MSE_C = ...;               % Take average over different channel realization for MSE
ASE_LMMSE_C = ...;         % Take average over different channel realization for ASE_LMMSE
ASE_LS_C = ...;            % Take average over different channel relaization for ASE_LS


figure
semilogy(...,...,'r-');
hold on
semilogy(...,...,'g--');
semilogy(...,...,'b-.');
hold off






