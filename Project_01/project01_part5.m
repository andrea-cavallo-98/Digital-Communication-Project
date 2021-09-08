%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #1
%
% Fifth step : use the companding technique with audio files
%

clear all
close all
clc

%% Read the audio file 
[Vin, Fs] = audioread('my_voice.ogg');

% Transpose Vin for the following calculations
Vin = Vin';

figure;

for Nbit = [4 6 8]

% Define constants
Nsamples = 1e6;
V = max(abs(Vin));
%Nbit = 6;
mu = 255;

M = 2^Nbit;
DeltaV = 2*V/M;

% Quantization of the signal
Partition = [-V+DeltaV:DeltaV:V-DeltaV]; % borders of intervals
CodeBook = [-V+DeltaV/2:DeltaV:+V-DeltaV/2]; % center elements of intervals

[Indexes, QuantizedSignal] = quantiz(Vin, Partition, CodeBook);

% Apply the compander (compressor with mu = 255)
Vin_compand = compand(Vin, mu, max(Vin), 'mu/compressor');
[Indexes_compand, QuantizedSignal_compand] = quantiz(Vin_compand, Partition, CodeBook);

% Convert decimal numbers to binary
bits_tx = de2bi(Indexes, Nbit);
bits_tx_compand = de2bi(Indexes_compand, Nbit);

%% Simulation

% The simulation is run for different values of P(e)
pe = logspace(-9, -1, 9);
for Counter=1:length(pe)

% Get bits at the receiver
bits_rx = bsc(bits_tx, pe(Counter));
bits_rx_compand = bsc(bits_tx_compand, pe(Counter));

% Convert binary to decimal and get output signals
IndexesOut = bi2de(bits_rx);
Vout = CodeBook(IndexesOut + 1);

% Apply the inverse companding component (expander with mu = 255)
IndexesOut_compand = bi2de(bits_rx_compand);
Vout_compand = CodeBook(IndexesOut_compand + 1);
Vout_compand = compand(Vout_compand, mu, max(Vout_compand), 'mu/expander');

% Overall error and SNR
e = Vout - Vin;
SNR_dB(Counter) = 10*log10(var(Vin)/var(e));

e_compand = Vout_compand - Vin;
SNR_dB_compand(Counter) = 10*log10(var(Vin)/var(e_compand));
end

h1 = plot(log10(pe), SNR_dB_compand, 'Marker', 'o');
set(h1, 'MarkerFaceColor', get(h1, 'color'));
hold on


end

xlabel('P(e)');
ylabel('SNR dB');
legend('4 bits', '6 bits', '8 bits')
title('Results with companding technique for different bits')


%%% Calculate theoretical SNR
pe_theory = logspace(-9, -1, 1e3);
SNR_dB_theory = 10*log10(M^2./(1+4*(M^2-1)*pe_theory));

%% Final plot

figure;
plot(log10(pe), SNR_dB, 'r', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r' );
hold on
plot(log10(pe), SNR_dB_compand, 'b', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_theory, 'k--');
xlabel('P(e)');
ylabel('SNR dB');
legend('Results with uniform quantization', 'Results with companding', 'Theoretical results for uniform quantization - uniform pdf')
title('Voice signal')
