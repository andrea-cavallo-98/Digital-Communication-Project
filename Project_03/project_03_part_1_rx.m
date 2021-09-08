%
% Andrea Cavallo
% matricola 245715
%
%
% PROJECT #3
% PART 1: Receiver
%
% NOTE: This script uses some variables from the tx file for error
% counting and filtering. Make sure to have those variables in the
% workspace when running this script.
%

clear all
close all
clc

%% General parameters

Fs = 44.1e3;
nBits = 8;
NumChannels = 1;
ID = -1;
l = 3;
f0 = 10e3;
record_loadn = 0; % 1 to record an audio file, 0 to used the stored sample
Tsim = 1/Fs;

%% Record audio data

if record_loadn
    recorder_obj = audiorecorder(Fs,nBits,NumChannels,ID);
    recordblocking(recorder_obj,l);
    y = getaudiodata(recorder_obj);
else
    y_load = load("tx_bits_64x64.mat"); 
    y = y_load.y;
    % This vector contains values from a recording I did and saved so that
    % I don't have to record every time
    % The vector "tx_bits.mat" contains the recording for the 10x10 image, the
    % vector "tx_bits_64x64.mat" contains the recording for the 64x64 image
end


%% Demodulation

s = ones(1, Ns);

t = 0:Tsim:Tsim*length(y)-Tsim;

% Multiply by cosine
y_dem = y .*cos(2*pi*f0*t)';

% Matched filter
H = (1./Ns).*conj(fftshift(fft(s,size(y_dem, 1))));
Y = fftshift(fft(y_dem));
R = H.*Y';
r = ifft(ifftshift(R));

%% Spectrum of received and filtered signal

[Freq_plot_r, PSDr] = project_03_myBartlett(r, round(size(r, 2)/10000), Fsim, 0);
[Freq_plot_ydem, PSDydem] = project_03_myBartlett(y_dem, round(size(y_dem, 1)/10000), Fsim, 0);
[Freq_plot_y, PSDy] = project_03_myBartlett(y, round(size(y, 1)/10000), Fsim, 0);
figure;
plot(Freq_plot_r / 1000, 10*log10(PSDr), 'r', 'LineWidth', 1);
hold on
plot(Freq_plot_y / 1000, 10*log10(PSDydem), 'k', 'LineWidth', 1);
hold on
plot(Freq_plot_y / 1000, 10*log10(PSDy), 'b', 'LineWidth', 1);
xlabel('Freq [kHz]');
ylabel('Power Spectral Density [dB]');
legend('Filtered signal', 'Demodulated signal', 'Received signal');
grid on
title('Spectra of signals at the receiver')

%% Synchronization

Barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
preamble = rectpulse(Barker,3);

% Normalization
y_dem_n = y_dem / max(abs(y_dem));

% Evaluate correlation
MFC_taps = rectpulse(flipud(preamble),SpS);
Corr = xcorr(MFC_taps, y_dem_n);

%% Plot correlation

figure;
plot(abs(Corr), 'k');
xlabel('Time');
ylabel('Amplitude');
title('Correlation with preamble');

%% Extract payload

% Find peaks
n1 = find(abs(Corr) == max(abs(Corr)));
n2_t = find(abs(Corr) == max([abs(Corr(1:n1-1000)) abs(Corr(n1+1000:end))]));

% Swap n1 and n2 if necessary
if n2_t < n1
    temp = n2_t;
    n2_t = n1;
    n1 = temp;
end    

% Exclude the second preamble 
n2 = n2_t - length(preamble) * SpS - 1;

% Extract the payload
payload = r(n1:n2);

% Normalization
payload_n = payload / max(abs(payload));

%% Eye diagram

lines = eyePlot(payload_n', 100, Ns, 25 * Ns,2);
grid on
xlabel('Time');
ylabel('Amplitude');
title('Eye diagram');
hold on
% Plot the average value
plot(linspace(-1,1,100), mean(mean(lines)) * ones(100, 1), 'k--')


%% Method 1: Calculate optimum sampling instant from eye diagram

lines_shifted = lines - mean(mean(lines)); % subtract average
maxdiff = 0;
for count = 1:size(lines,2)

    % Minimum positive value and maximum negative value of eye
    min_pos = min(lines_shifted(lines_shifted(:,count) > 0,count));
    max_neg = max(lines_shifted(lines_shifted(:,count) < 0,count));
    diff = abs(min_pos - max_neg);
    % Identify the sampling instant of maximum eye opening
    if diff > maxdiff
        kopt = count;
        maxdiff = diff;
    end
    
end
kopt

%% Method 1: Sampling and decision

min_err = length(MessageBits);

% Threshold voltage
Vth = abs(mean(mean(lines))); % Average value of eye diagram as threshold

% Decision at optimum sampling instant
r_samp = payload_n(kopt:Ns:end);
r_final(abs(r_samp) > Vth) = 1;
r_final(abs(r_samp) <= Vth) = 0;

% BER counting using original vector
err = sum(flip(r_final) ~= MessageBits');
BER = err / length(MessageBits)
    
%% Scatter plot

plot(r_samp_k(flip(MessageBits) == 0), zeros(size(r_samp(flip(MessageBits) == 0))), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'LineStyle', 'None');
hold on
plot(r_samp_k(flip(MessageBits) == 1), zeros(size(r_samp(flip(MessageBits) == 1))), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'LineStyle', 'None');
hold on
% Add a circle whose radius is the selected threshold
%t_var = [0:0.1:6.5]; 
%plot(Vth * cos(t_var), Vth * sin(t_var), 'k', 'LineWidth', 1)
%plot(-0.3,0, 'o');
grid on
xlabel('Real part');
ylabel('Imaginary part');
legend('Samples when sent bit is 0', 'Samples when sent bit is 1', 'Threshold');
title('Scatter plot of samples at optimum sampling time');

%% Method 2: K-Means clustering
% The Statistics and Machine Learning Toolbox is required 

max_min_dist = 0;

for k = 1:Ns % test on all possible sampling instants
    r_samp_k = payload_n(k:Ns:end); % Extract samples at given sampling instant
    [idx_k, C, sumd] = kmeans(r_samp_k', 2); % Apply clustering algorithm
    
    % The minimum distance between the clusters is the distance between two
    % extreme point of each cluster (for this particular 1D configuration)
    min_dist = abs(min(r_samp_k(idx_k == 1)) - max(r_samp_k(idx_k == 2)));
    if abs(max(r_samp_k(idx_k == 1)) - min(r_samp_k(idx_k == 2))) < min_dist
        min_dist = abs(max(r_samp_k(idx_k == 1)) - min(r_samp_k(idx_k == 2)));
    end
    
    % The optimum sampling time is chosen according to the maximum minimum
    % distance between the clusters
    if min_dist > max_min_dist
       max_min_dist = min_dist;
       idx_k_best = idx_k;
       C_best = C;
       kopt = k;
       r_samp_k_best = r_samp_k;
    end
    
end

% The centroids are compared to assign the transmitted value to each
% cluster
if abs(C_best(1)) > abs(C_best(2))
    r_final(idx_k_best == 1) = 1;
    r_final(idx_k_best == 2) = 0;
else
    r_final(idx_k_best == 1) = 0;
    r_final(idx_k_best == 2) = 1;
end

% Errors and BER counting
err = sum(flip(r_final) ~= MessageBits');
BER = err / length(MessageBits)

%% Method 2: K-Means clustering plot

plot(r_samp_k_best(idx_k_best == 1), zeros(size(r_samp_k_best(idx_k_best == 1))), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'LineStyle', 'None');
hold on
plot(r_samp_k_best(idx_k_best == 2), zeros(size(r_samp_k_best(idx_k_best == 2))), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'LineStyle', 'None');
hold on
plot(C_best, zeros(size(C_best)), 'Marker', 'x', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineStyle', 'None', 'LineWidth', 2)
grid on
xlabel('Real part');
ylabel('Imaginary part');
legend('Cluster 1', 'Cluster 2', 'Centroids');
title('Clustering diagram of samples at optimum sampling time');

%% Rebuild image

Abin_r = reshape(uint8(flip(r_final)), Size, 8);
A_r = bi2de(Abin_r, 'left-msb');
DataRxImg = reshape(A_r, Rows, Columns, D);
figure;
image(DataRxImg);
title('Image at receiver');



