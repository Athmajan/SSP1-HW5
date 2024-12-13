clear all; close all; clc;

MC = 1000;                      % number of Monte Carlo runs for data
MC_channel = 1000;              % number of Monte Carlo runs for channel
N_TX = 2;                       % number of TX antennas
N_RX = 2;                       % number of RX antennas

EbN0dB = 0:0.1:10;              % SNR vector in dB
NoSNRPoints = length(EbN0dB);   % length of SNR vector
EbN0 = 10.^(EbN0dB./10);        % SNR in linear scale
C = 1./(2*EbN0);                % Find noise variance per I and Q component
sig = sqrt(C) ;                 % Find noise standard deviation per I and Q component

A = sign(randn(N_TX,MC));       % random data for all runs, the same for all SNR points
C_A = diag(var(A'));                        % Find data variance (variance of A)
ASE_LS = zeros(NoSNRPoints,MC_channel);     % initialize the LS average square error (ASE)
ASE_LMMSE = zeros(NoSNRPoints,MC_channel);  % LMMSEE ASE
MSE_LMMSE = zeros(NoSNRPoints,MC_channel);  % theoretical MSE for LMMSE

LS_BER = zeros(NoSNRPoints,1);    % LS average square error (ASE)
LMMSE_BER = zeros(NoSNRPoints,1); % LMMSEE ASE

for i = 1:NoSNRPoints        % generate noise realization for SNR values using
    for c = 1:MC_channel     %generate channel realization 
        H = randn(N_RX, N_TX) + 1i*randn(N_RX, N_TX);              % generate complex channel response with total variance 1
        w = sig(i).*(randn(N_TX,MC) + 1i*randn(N_TX,MC));          % generate complex noise with different sig for different SNR        
        X = H*A+ w;
        
        C_w = eye(N_TX)*C(i);      
        
        LS = inv(H'*inv(C_w)*H)*H'*inv(C_w)*X;               % Implement LS estimator and take only real part which contain the data
        LMMSE = inv(inv(C_A) + H'*inv(C_w)*H)*H'*inv(C_w)*X; % Implement LMMSE estimator and take only real part which contain the data
        ASE_LS(i,c) = mean((abs(LS-A)).^2,'all');            % Calculate average square error for LS
        ASE_LMMSE(i,c) = mean((abs(LMMSE-A)).^2,'all');                         % Calculate average square error for LMMSE
        MSE_LMMSE(i,c) = sum(abs(diag(real(inv(inv(C_A) + H'*inv(C_w)*H)))))/2; % Implement theoritical MSE for LMMSE by taking real parts of diagonal elements
        
        LS_BER(i) = LS_BER(i)+ sum((sign(real(LS))~=A),'all');
        LMMSE_BER(i) = LMMSE_BER(i)+ sum((sign(real(LMMSE))~=A),'all');    
    end
end

MSE_C = mean(MSE_LMMSE,2);            % Take average over different channel realization for MSE
ASE_LMMSE_C = mean(ASE_LMMSE,2);      % Take average over different channel realization for ASE_LMMSE
ASE_LS_C = mean(ASE_LS,2);            % Take average over different channel relaization for ASE_LS

LS_BER = LS_BER/(N_TX*MC_channel*MC);
LMMSE_BER = LMMSE_BER/(N_TX*MC_channel*MC);

f1 = figure
plot(EbN0dB,MSE_C,'b-');
hold on
plot(EbN0dB,LS_BER,'g-');
plot(EbN0dB,LMMSE_BER,'r-');
hold off
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Theortical MSE', 'BER of LSE', 'BER of LMMSE');

f2 = figure
plot(EbN0dB,MSE_C,'b-');
hold on
plot(EbN0dB,ASE_LMMSE_C,'g-');
hold off
xlabel('SNR (dB)');
ylabel('Average Absolute Square Error (ASE SE)');
legend( 'Theoratical MSE', 'LMMSE');
%

