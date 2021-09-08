%
% Andrea Cavallo
% matricola 245715
%
%
% PROJECT #3
% PART 1: Transmitter
%


clear all
close all
clc

%% Create image

image_10x10 = 0; % 1 to use the 10x10 image, 0 to use the 64x64 image

% I created a small image (10x10) to make transmission faster
if image_10x10 
    Mat = [0 10 20 30 40 50 60 70 80 90; 5 15 25 35 45 55 65 75 85 95;
        0 10 20 30 40 50 60 70 80 90; 5 15 25 35 45 55 65 75 85 95;
        0 10 20 30 40 50 60 70 80 90; 5 15 25 35 45 55 65 75 85 95;
        0 10 20 30 40 50 60 70 80 90; 5 15 25 35 45 55 65 75 85 95;
        0 10 20 30 40 50 60 70 80 90; 5 15 25 35 45 55 65 75 85 95];
else
    Mat = uint8(get(0,'DefaultImageCData')); % Default Matlab image, 64x64
end

%% Read and convert image

DataTxImg = Mat;
figure;
image(Mat);
title('Image at transmitter');
[Rows,Columns,D] = size(DataTxImg );
A = reshape(DataTxImg,Rows*Columns*D,1); 
Abin = de2bi(A,8,'left-msb'); 
Size = size(Abin,1); 
MessageBits = reshape(Abin,Size*8,1);

%% Add preamble

Barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
preamble = rectpulse(Barker,3);
dummy = zeros(1, 10000);
frame = [dummy preamble double(MessageBits') preamble dummy];

%% General parameters

Fsim = 44.1e3;
f0 = 10e3;
M = 2;
Ns = 8;
SpS = Ns * log2(M);
Rs = Fsim / SpS;
Rb = Rs * log2(M);
Nbits = length(frame);
Nsamples = Nbits*Ns;
Tsim = 1/Fsim;
Fsound = 20e3;

%% Modulated signal (2PAM)

s = ones(1, Ns);

x_i = zeros(Nsamples,1);
for n = 1:Nbits
   idx = (n-1)*Ns +(1:Ns);
   x_i(idx,1) = (frame(n).*s)'; % unipolar configuration
end

t = [0:length(x_i)-1]*Tsim;

x = x_i.*cos(2*pi*f0*t)';

%% Spectrum of transmitted signal

[Freq_plot_i, PSDx_i] = project_03_myBartlett(x_i, round(size(x_i, 1)/10000), Fsim, 0);
[Freq_plot, PSDx] = project_03_myBartlett(x, round(size(x, 1)/10000), Fsim, 0);
figure;
plot(Freq_plot, 10*log10(PSDx), 'r', 'LineWidth', 1);
hold on
plot(Freq_plot_i, 10*log10(PSDx_i), 'k', 'LineWidth', 1);
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB]');
legend('Modulated signal', 'Baseband signal');
grid on
title('Spectrum of baseband and modulated signal at the transmitter');

%% Generate audio file

% Normalization
x_n = x / abs(max(x));
audiowrite('my_imagefast.wav', x_n, Fsim);


