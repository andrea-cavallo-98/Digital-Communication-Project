%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #2
% 2-PAM system
%
% DESCRIPTION:
% This file contains the simulation of 2-PAM system for several filters
% and pulse shapes. It also computes the Eb/No penalty when using
% an RC filter.
%


function [Freq_plot, PSDx, EbNo_dB, BERth, BER] = project_02_2PAM( tx_filter_type, rx_filter_type, Nbits, Rb, Ns, roll_off, f3dB_coeff )


%% General parameters 

Tb = 1./Rb;                      % Time Bit [s]
Nsamples = Ns*Nbits;             % Number of samples
Fsim = Ns*Rb;                    % Simulation Bandwidth [Hz]
Tsim = 1./Fsim;                  % Sample time [s]
EbNo_dB = [2:10];                % Eb/No [dB]
EbNo = 10.^(EbNo_dB*0.1);        % Eb/No
BERth = 0.5*erfc(sqrt(EbNo));    % Theoretical BER
Rs = Rb;                         % Symbol Rate
f3dB = f3dB_coeff * Rs;          % Bandwidth for the RC filter
nfft = Nsamples;                 % Samples for Fourier transform
stepFreq = Fsim/nfft;            % Step-frequency
maxFreq = +Fsim/2-stepFreq;      % Max frequency
minFreq = -Fsim/2;               % Min frequency
Freq = (minFreq:stepFreq:maxFreq)'; % Frequency vector

%% Transmitter

% Bit generation 
Bits = randi([0 1], Nbits, 1);
% Antipodal representation
values(Bits == 0) = -1;
values(Bits == 1) = 1;

[x,H] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nbits, Fsim, Rs, f3dB, Freq, roll_off);

%% Spectral Analysis

% Spectrum of the transmitter signal

[Freq_plot, PSDx] = myBartlett(x', 500, Fsim);



%% Transmission and filtering

for ii = 1:length(EbNo_dB)
    %% AWGN channel
    Ps = var(x);            % Signal Power
    No = (Ps/Rb)./EbNo(ii);
    Pn = No/2*Fsim;
    noise = sqrt(Pn).*randn(Nsamples,1); % noise signal
    y = x+noise; % add the WGN to the signal
    
    %% Apply filter
    
    % In frequency domain
    Y = fftshift(fft(y));
    R = H.*Y;
    r = ifft(ifftshift(R));
    
    %% Ber counting and optimum sampling time
    
    % Threshold voltage
    Vth = 0.0;

    for k = 1:Ns
        rk = r(k:Ns:end);
        rkd = (rk>Vth);
        err(k) =sum(rkd~=Bits);
    end
    
    [BER(ii),kopt]= min(err./Nbits);
    
end
end