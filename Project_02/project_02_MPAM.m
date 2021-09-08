%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #2
% M-PAM system
%
% NOTE:
% This file contains the simulation of M-PAM system (4 or 8) for several filters
% and pulse shapes. It also computes the Eb/No penalty when using an RC
% filter.
%


function [Freq_plot, PSDx, EbNo_dB, theoryBer, theorySer, BER, SER] = project_02_MPAM(M, GrayCoding, tx_filter_type, rx_filter_type, roll_off, Nsymbols, Ns, Rb, f3dB_coeff)

%% General parameters

Rs = Rb / log2(M);                                              % symbol rate
Nbits = Nsymbols * log2(M);                                     % number of bits
[values_range, bit_mapping] = create_codings(M, GrayCoding);    % values and bit mapping
EbNo_dB = 4+M:10+M;                                              % Eb/No [db]
Nsamples = Ns*Nsymbols;                                         % Number of samples
Fsim = Ns*Rs;                                                   % Simulation Bandwidth [Hz]
EbNo = 10.^(EbNo_dB*0.1);                                       % Eb/No      
nfft = Nsamples;                                                % Samples for Fourier transform 
stepFreq = Fsim/nfft;                                           % Step-frequency
maxFreq = +Fsim/2-stepFreq;                                     % Max frequency
minFreq = -Fsim/2;                                              % Min frequency
Freq = (minFreq:stepFreq:maxFreq)';                             % Frequency vector
f3dB = f3dB_coeff * Rs;                                         % Bandwidth for RC filter

%% Transmitter

% Generate values amongst values_range
values = randsrc(1,Nsymbols,values_range);

% Generate string of bits corresponding to the vector values_range
values_bits = bit_mapping(getindex(values(1),M), :);
for count = 2:size(values,2)
   values_bits = [values_bits bit_mapping(getindex(values(count),M), :)]; 
end 
values_bits = values_bits';

[x,H] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nsymbols, Fsim, Rs, f3dB, Freq, roll_off);

%% Spectral Analysis 

% Spectrum of transmitter signal
[Freq_plot, PSDx] = myBartlett(x', 500, Fsim);


%% Transmission and filtering

for ii = 1:length(EbNo_dB)
    %% AWGN channel
    Ps = var(x);            % Signal Power
    No = (Ps/Rb)./EbNo(ii);
    Pn = No/2*Fsim;
    noise = sqrt(Pn).*randn(Nsamples,1); % noise signal
    y = x+noise; % add WGN to transmitter signal

    %% Apply filter
    
    % In frequency domain
    Y = fftshift(fft(y));
    R = H.*Y;
    r = ifft(ifftshift(R));
    
    %% Ser and Ber counting and optimum sampling time

    % SER
    min_err = 1e10; % needed for BER
    
    for k = 1:Ns
        rk = r(k:Ns:end);
        rkd = decode_values(rk,M); % get symbols from signal amplitude
        err(k) = sum(rkd~=values);
        if err(k) < min_err
           min_err = err(k);
           rkd_opt = rkd;
        end
    end
    
    [SER(ii),kopt]= min(err./Nsymbols);
    
    % BER
    errors = 0;
    c = 0;
    for count = 1:size(rkd_opt,2)
       if rkd_opt(count) ~= values(count)
           % if the symbol is wrong, the Hamming distance between the two
           % symbols is computed
           errors = errors + sum(bit_mapping(getindex(rkd_opt(count),M),:) ~= bit_mapping(getindex(values(count),M),:));
       end
    end
    
    BER(ii)  = errors / Nbits;
    
end

%% Find theoretical SER and BER

% SER
m = log2(M);
g = 3*m/(M^2-1);
theorySer = (M-1)/M*erfc(sqrt(g*EbNo));

% BER
BER1 = 0;
for j=0:M-2
    dH = sum(xor(bit_mapping(j+2,:),bit_mapping(j+1,:)));
    BER1 = (1/m/M)*dH.*erfc(sqrt(g*EbNo)) + BER1;
end

BER2 = 0;
for j=0:M-3
    for k = j+1:M-2
        dH1 = sum(xor(bit_mapping(k+2,:),bit_mapping(j+1,:)));
        dH2 = sum(xor(bit_mapping(k+1,:),bit_mapping(j+1,:)));
        BER2 = (1/m/M)*(dH1-dH2).*erfc(sqrt((2*k+2-2*j-2+1)^2*g*EbNo)) + BER2;
    end
end

theoryBer = BER1 + BER2; 


end

%% Functions

% Create the possible values and the bit mapping for different M values
function [values_range, bit_mapping] = create_codings(M, GrayCoding)

if M == 4
    values_range = [-3 -1 1 3];
    if GrayCoding
        bit_mapping = [0 0; 0 1; 1 1; 1 0];
    end
    if ~GrayCoding
        bit_mapping = [0 0; 0 1; 1 0; 1 1];
    end
end

if M == 8
    values_range = [-7 -5 -3 -1 1 3 5 7];
    if GrayCoding
        bit_mapping = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1; 1 1 1];
    end
    if ~GrayCoding
        bit_mapping = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
    end
end

end


% Get symbols from signal amplitude
function rkd = decode_values(rk, M)

if M == 4
    rkd(rk< -2) = -3;
    rkd(rk>= 2) = 3;
    rkd(rk>=-2 & rk<0) = -1;
    rkd(rk>=0 & rk<2) = 1;
end

if M == 8
    rkd(rk< -6) = -7;
    rkd(rk>= 6) = 7;
    rkd(rk< -4 & rk >= -6) = -5;
    rkd(rk>= 4 & rk < 6) = 5;
    rkd(rk< -2 & rk >= -4) = -3;
    rkd(rk>= 2 & rk < 4) = 3;
    rkd(rk>=-2 & rk<0) = -1;
    rkd(rk>=0 & rk<2) = 1;
end

end

% Get index in the bit mapping table from the symbol
function idx = getindex(a,M)
    if M == 4
        idx = (a + 5) / 2;
    end
    if M == 8
        idx = (a + 9) / 2;
    end
end

