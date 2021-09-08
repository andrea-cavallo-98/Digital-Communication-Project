%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #1
%
% First step : test the simulator with different values of bits
% transmitted and different probability of error in the channel
%


clear all
close all
clc

Nsamples = 1e6;
V = 8;
Nbit_range = [8 6 4];

% creating some containers to store data and plot at the end
SNR_dB_theory = zeros(3, 1e3);
SNR_dB_vector = zeros(3,9);
pe_vector = zeros(3,9);
SNR_dB_bsceffect_vector = zeros(3,9);

for Nbit = Nbit_range
M = 2^Nbit;
DeltaV = 2*V/M;

% generate the sample with rand function
Vin = 2*V*rand(1, Nsamples) - V;

% Quantization of the signal
Partition = [-V+DeltaV:DeltaV:V-DeltaV]; % borders of intervals
CodeBook = [-V+DeltaV/2:DeltaV:+V-DeltaV/2]; % center elements of intervals

[Indexes, QuantizedSignal] = quantiz(Vin, Partition, CodeBook);

% convert decimal numbers to binary
bits_tx = de2bi(Indexes, Nbit);

% The simulation is run for different values of P(e)
pe = logspace(-9, -1, 9);
for Counter=1:length(pe)

% get bits at the receiver
bits_rx = bsc(bits_tx, pe(Counter));

% Convert to decimal and get Vout
IndexesOut = bi2de(bits_rx);
Vout = CodeBook(IndexesOut + 1);

%%% error and SNR due to BSC
eB = Vout - QuantizedSignal;
SNR_dB_bsceffect(Counter)=  10*log10(var(Vin)/var(eB));

% total error and SNR
e = Vout - Vin;
SNR_dB(Counter) = 10*log10(var(Vin)/var(e));
end

%%% calculate theoretical curve with good precision
pe_theory = logspace(-9, -1, 1e3);
SNR_dB_theory(Nbit/2-1, :) = 10*log10(M^2./(1+4*(M^2-1)*pe_theory));
SNR_dB_vector(Nbit/2-1, :) = SNR_dB;
pe_vector(Nbit/2-1, :) = pe;
SNR_dB_bsceffect_vector(Nbit/2-1,:) = SNR_dB_bsceffect;
end

%%% Plots of the signal
% plot the signal and its pdf
figure;
plot(Vin(1:50));
title('Transmitted signal');
xlabel('Time');
ylabel('Vin');

figure;
histogram(Vin, 200, 'Normalization', 'pdf')
title('PDF of transmitted signal');
xlabel('x');
ylabel('PDF');


%% Spectrum

[f_Bartlett, F_Bartlett] = myBartlett(Vin, Nsamples/1000, Nsamples);
figure;
plot(f_Bartlett, 20 * log10(F_Bartlett));
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB]');
title('PSD of transmitted signal');

%% Plot to visualize the quantization 

figure;
plot(Vin(1:50),'r-');
hold on;
grid on;
plot(QuantizedSignal(1:50), 'bo-');
xlabel('Time');
ylabel('Signals [V]');
legend('Original signal', 'Quantized signal');
title('Quantization');

% Compare the original signal to the signal at the receiver
figure;
plot(Vin(1:50), 'r');
hold on;
grid on;
plot(Vout(1:50), 'g');
xlabel('Time');
ylabel('Signals [V]');
legend('Transmitted signal', 'Received signal');
title('Compare transmitted and received signal');

%% Final plots

%%% Final plot for different N_bits
figure;
plot(log10(pe_vector(1,:)), SNR_dB_vector(1,:), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'None');
hold on
plot(log10(pe_vector(2,:)), SNR_dB_vector(2,:), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineStyle', 'None');
hold on
plot(log10(pe_vector(3,:)), SNR_dB_vector(3,:), 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineStyle', 'None');
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_theory(1,:), 'r', log10(pe_theory), SNR_dB_theory(2,:), 'b', log10(pe_theory), SNR_dB_theory(3,:), 'k');
hold on;
xlabel('P(e)');
ylabel('SNR dB');
legend('N bit = 4', 'N bit = 6', 'N bit = 8', 'N bit = 4', 'N bit = 6', 'N bit = 8');

%%% Plot with BSC error for N_bit = 6
figure;
plot(log10(pe_vector(2,:)), SNR_dB_vector(2,:), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'None' );
hold on
plot(log10(pe_vector(2,:)), SNR_dB_bsceffect_vector(2,:), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineStyle', 'None');
grid on;
hold on;
SNR_dB_bsc_theory = 10*log10(1./(4*pe_theory));
plot(log10(pe_theory), SNR_dB_theory(2,:), 'r', log10(pe_theory), SNR_dB_bsc_theory, 'b')
xlabel('P(e)');
ylabel('SNR dB');
legend('Experimental SNR', 'Experimental SNR due to BSC effect', 'Theoretical SNR', 'Theoretical SNR due to BSC effect')

