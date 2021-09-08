%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #1
%
% Third step : perform transmission analysis from two audio files 
% (music and human voice)
%

clear all
close all
clc

%% Read files and plot spectra and PDF

% Read the file with human voice
[my_voice, Fs_my_voice] = audioread('my_voice.ogg');
my_voice = my_voice'; % the vector needs to be transposed for next operations

%%% Plot the pdf of human voice
figure;
histogram(my_voice, 200, 'Normalization', 'pdf')
xlabel('x');
ylabel('PDF')
title('PDF of human voice');

figure;
[f_Bartlett, F_Bartlett] = myBartlett(my_voice, round(length(my_voice)/1000), length(my_voice));
plot(f_Bartlett, 20 * log10(F_Bartlett));
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB]');
title('PSD of human voice');

% Read the file with music
[music, Fs_music] = audioread('we_are_the_champions.ogg');
music = music'; % the vector needs to be transposed for next operations

%%% Plot the pdf of music
figure;
histogram(music, 200, 'Normalization', 'pdf')
xlabel('x');
ylabel('PDF')
title('PDF of music');

figure;
[f_Bartlett, F_Bartlett] = myBartlett(music, round(length(music)/1000), length(music));
plot(f_Bartlett, 20 * log10(F_Bartlett));
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB]');
title('PSD of music');

%% Simulation

% Define constants for the simulation
Nsamples = 1e6;
V_my_voice = max(abs(my_voice)); % max value of the sample of voice (avoid clipping)
V_music = max(abs(music)); % max value of the sample of music (avoid clipping)
Nbit = 6;

M = 2^Nbit;
DeltaV_my_voice = 2*V_my_voice/M;
DeltaV_music = 2*V_music/M;

% Quantization of the samples
Partition_my_voice = [-V_my_voice+DeltaV_my_voice:DeltaV_my_voice:V_my_voice-DeltaV_my_voice]; % borders of intervals
CodeBook_my_voice = [-V_my_voice+DeltaV_my_voice/2:DeltaV_my_voice:+V_my_voice-DeltaV_my_voice/2]; % center elements of intervals

Partition_music = [-V_music+DeltaV_music:DeltaV_music:V_music-DeltaV_music]; % borders of intervals
CodeBook_music = [-V_music+DeltaV_music/2:DeltaV_music:+V_music-DeltaV_music/2]; % center elements of intervals

[Indexes_my_voice, QuantizedSignal_my_voice] = quantiz(my_voice, Partition_my_voice, CodeBook_my_voice);
[Indexes_music, QuantizedSignal_music] = quantiz(music, Partition_music, CodeBook_music);

% Convert decimal numbers to binary
bits_tx_my_voice = de2bi(Indexes_my_voice, Nbit);
bits_tx_music = de2bi(Indexes_music, Nbit);

% The simulation is run for different values of P(e)
pe = logspace(-9, -1, 9);
for Counter=1:length(pe)

% Get bits at the receiver
bits_rx_my_voice = bsc(bits_tx_my_voice, pe(Counter));
bits_rx_music = bsc(bits_tx_music, pe(Counter));

% Convert binary to decimal and get output signals
IndexesOut_my_voice = bi2de(bits_rx_my_voice);
Vout_my_voice = CodeBook_my_voice(IndexesOut_my_voice + 1);

IndexesOut_music = bi2de(bits_rx_music);
Vout_music = CodeBook_music(IndexesOut_music + 1);

% Overall errors and SNRs
e_my_voice = Vout_my_voice - my_voice;
SNR_dB_my_voice(Counter) = 10*log10(var(my_voice)/var(e_my_voice));

e_music = Vout_music - music;
SNR_dB_music(Counter) = 10*log10(var(music)/var(e_music));
end

%%% Calculate theoretical SNR 
pe_theory = logspace(-9, -1, 1e3);
SNR_dB_theory = 10*log10(M^2./(1+4*(M^2-1)*pe_theory));

%%% Estimation of the loss coefficients 
Vin_my_voice = 2*V_my_voice*rand(1,Nsamples) - V_my_voice; % random sample with uniform pdf
a_my_voice = 10*log10(var(Vin_my_voice)./var(my_voice)); 

Vin_music = 2*V_music*rand(1,Nsamples) - V_music; % random sample with uniform pdf
a_music = 10*log10(var(Vin_music)./var(music));

%% Final plot

figure;
plot(log10(pe), SNR_dB_my_voice, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'None');
hold on
plot(log10(pe), SNR_dB_music, 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineStyle', 'None')
grid on;
hold on;
plot(log10(pe_theory), SNR_dB_theory, 'k--', log10(pe_theory), SNR_dB_theory - a_my_voice,'r', log10(pe_theory), SNR_dB_theory - a_music, 'b');
hold on;
xlabel('P(e)');
ylabel('SNR dB');
legend('Human voice signal', 'Music signal', 'Theoretical value for signal with uniform pdf', 'Theoretical value for voice signal', 'Theoretical value for music signal');

