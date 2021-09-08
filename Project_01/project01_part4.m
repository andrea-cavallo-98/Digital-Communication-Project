%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #1
%
% Fourth step : use non uniform quantization with audio files
%

clear all
close all
clc

%% Read the audio file 
[Vin, Fs] = audioread('my_voice.ogg');

% Define constants
V = max(abs(Vin)); % max of the sample (avoid clipping)
Nbit = 6;

M = 2^Nbit;
DeltaV = 2*V/M;

% Transpose Vin to perform the following calculations
Vin = Vin';

% Quantization of the signal using both uniform and non uniform quantization
% (non uniform quantization is achieved using Lloyd's algorithm)
Partition = [-V+DeltaV:DeltaV:V-DeltaV]; % borders of intervals
CodeBook = [-V+DeltaV/2:DeltaV:+V-DeltaV/2]; % center elements of intervals

[Partition_improved,CodeBook_improved] = lloyds(Vin,CodeBook);

[Indexes, QuantizedSignal] = quantiz(Vin, Partition, CodeBook);
[Indexes_improved, QuantizedSignal_improved] = quantiz(Vin, Partition_improved, CodeBook_improved);

% Convert decimal numbers to binary
bits_tx = de2bi(Indexes, Nbit);
bits_tx_improved = de2bi(Indexes_improved, Nbit);

%% Simulation

% The simulation is run for different values of P(e)
pe = logspace(-9, -1, 9);
for Counter=1:length(pe)

% Get bits at the receiver
bits_rx = bsc(bits_tx, pe(Counter));
bits_rx_improved = bsc(bits_tx_improved, pe(Counter));

% Convert binary to decimal and get outputs
IndexesOut = bi2de(bits_rx);
Vout = CodeBook(IndexesOut + 1);

IndexesOut_improved = bi2de(bits_rx_improved);
Vout_improved = CodeBook_improved(IndexesOut_improved + 1);

% Overall errors and SRNs
e = Vout - Vin;
SNR_dB(Counter) = 10*log10(var(Vin)/var(e));

e_improved = Vout_improved - Vin;
SNR_dB_improved(Counter) = 10*log10(var(Vin)/var(e_improved));
end

%%% Calculate theoretical SNR
pe_theory = logspace(-9, -1, 1e3);
SNR_dB_theory = 10*log10(M^2./(1+4*(M^2-1)*pe_theory));

%% Final plot

figure;
plot(log10(pe), SNR_dB, 'r', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
hold on
plot(log10(pe), SNR_dB_improved, 'b', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_theory, 'k--');
xlabel('P(e)');
ylabel('SNR dB');
legend('Results with uniform quantization', 'Results with non uniform quantization', 'Results with uniform quantization for signal with uniform pdf')
title('Voice signal');
