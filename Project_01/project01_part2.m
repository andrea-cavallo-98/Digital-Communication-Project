%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #1
%
% Second step : check with ADCs that use as output values the lower 
% and the upper end of the quantization intervals
%

clear all
close all
clc

Nsamples = 1e6;
V = 8;
Nbit = 6; % the simulation is run for a channel with 6 bits

M = 2^Nbit;
DeltaV = 2*V/M;

% generate the sample with rand function
Vin = 2*V*rand(1, Nsamples) - V;

% Quantization of the signal with different CodeBooks
Partition = [-V+DeltaV:DeltaV:V-DeltaV]; % borders of intervals
CodeBook_lower = [-V:DeltaV:+V-DeltaV]; % lower elements of intervals
CodeBook_upper = [-V+DeltaV:DeltaV:+V]; % upper elements of intervals

[Indexes_lower, QuantizedSignal_lower] = quantiz(Vin, Partition, CodeBook_lower);
[Indexes_upper, QuantizedSignal_upper] = quantiz(Vin, Partition, CodeBook_upper);

%%% Plot the pdfs of quantized signal (as expected they are translated)
figure;
histogram(QuantizedSignal_lower, 64, 'Normalization', 'pdf');
xlabel('x');
ylabel('PDF');
title('PDF of signal quantized with lower elements of intervals')
figure;
histogram(QuantizedSignal_upper, 64, 'Normalization', 'pdf');
xlabel('x');
ylabel('PDF');
title('PDF of signal quantized with lower elements of intervals')

% converts decimal numbers to binary
bits_tx_lower = de2bi(Indexes_lower, Nbit);
bits_tx_upper = de2bi(Indexes_upper, Nbit);

% The simulation is run for different values of P(e)
pe = logspace(-9, -1, 9);
for Counter=1:length(pe)

% get bits at the receiver
bits_rx_lower = bsc(bits_tx_lower, pe(Counter));
bits_rx_upper = bsc(bits_tx_upper, pe(Counter));

% Convert to decimal and get Vout
IndexesOut_lower = bi2de(bits_rx_lower);
IndexesOut_upper = bi2de(bits_rx_upper);
Vout_lower = CodeBook_lower(IndexesOut_lower + 1);
Vout_upper = CodeBook_upper(IndexesOut_upper + 1);

%%% Errors and SNRs due to BSC
eB_lower = Vout_lower - QuantizedSignal_lower;
eB_upper = Vout_upper - QuantizedSignal_upper;
SNR_dB_bsceffect_lower(Counter) =  10*log10(var(Vin)/var(eB_lower));
SNR_dB_bsceffect_upper(Counter) =  10*log10(var(Vin)/var(eB_upper));

% Overall errors and SNRs
e_lower = Vout_lower - Vin;
e_upper = Vout_upper - Vin;

% Pn is not the variance!!!! Because the signal is not zero average
% We need to use the mean of the error
SNR_dB_lower(Counter) = 10*log10(var(Vin)/mean(e_lower.^2));
SNR_dB_upper(Counter) = 10*log10(var(Vin)/mean(e_upper.^2));
end

% max eq = DeltaV
% SNR = M^2/4

%%% Calculate theoretical SNR with good precision
pe_theory = logspace(-9, -1, 1e3);
SNR_dB_theory_uniform_quant = 10*log10(M^2./(1+4*(M^2-1)*pe_theory));
SNR_dB_bsc_theory = 10*log10(1./(4*pe_theory));
% The theoretical SNR changes: M^2/4
SNR_dB_quant_theory = 10*log10(M^2./4);


%% Final plot

figure;
plot(log10(pe), SNR_dB_lower,  'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'None');
hold on
plot(log10(pe), SNR_dB_bsceffect_lower,  'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineStyle', 'None')
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_quant_theory*ones(size(pe_theory)), 'r', log10(pe_theory), SNR_dB_bsc_theory, 'k', log10(pe_theory), SNR_dB_theory_uniform_quant, 'b');
xlabel('P(e)');
ylabel('SNR dB');
legend('SNR with lower elements of intervals',...
    'SNR due to BSC with lower elements of intervals',...
     'SNR due to quantization error', 'Expected SNR due to BSC', ...
    'SNR for center elements of intervals');
title('Results for lower elements of intervals')

figure;
plot(log10(pe), SNR_dB_upper,  'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'None')
hold on
plot(log10(pe), SNR_dB_bsceffect_upper,  'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineStyle', 'None');
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_quant_theory*ones(size(pe_theory)), 'r', log10(pe_theory), SNR_dB_bsc_theory, 'k', log10(pe_theory), SNR_dB_theory_uniform_quant, 'b');
xlabel('P(e)');
ylabel('SNR dB');
legend('SNR with upper elements of intervals',...
    'SNR due to BSC with upper elements of intervals', 'SNR due to quantization error', 'Expected SNR due to BSC', ...
    'SNR for center elements of intervals');
title('Results for upper elements of intervals')
