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


function [penalty, f3dB_coef] = project_02_MPAM_RCpenalty(M, GrayCoding, tx_filter_type, roll_off, Nsymbols, Ns, Rb)

%% General parameters

Rs = Rb / log2(M);                                              % symbol rate
Nbits = Nsymbols * log2(M);                                     % number of bits
[values_range, bit_mapping] = create_codings(M, GrayCoding);    % values and bit mapping

% EbNo values are selected according to expected results
if M == 4 
    if tx_filter_type == "NRZ"
        EbNo_dB = 10:17;
    else
        EbNo_dB = 10:24;
    end
else
    if tx_filter_type == "NRZ"
        EbNo_dB = 10:27;
    else
        EbNo_dB = 10:26;
    end
end

Nsamples = Ns*Nsymbols;                                         % Number of samples
Fsim = Ns*Rs;                                                   % Simulation Bandwidth [Hz]
EbNo = 10.^(EbNo_dB*0.1);                                       % Eb/No      
nfft = Nsamples;                                                % Samples for Fourier transform 
stepFreq = Fsim/nfft;                                           % Step-frequency
maxFreq = +Fsim/2-stepFreq;                                     % Max frequency
minFreq = -Fsim/2;                                              % Min frequency
Freq = (minFreq:stepFreq:maxFreq)';                             % Frequency vector


%% Transmitter

% Generate values amongst values_range
values = randsrc(1,Nsymbols,values_range);

% Generate string of bits corresponding to the vector values_range
values_bits = bit_mapping(getindex(values(1),M), :);
for count = 2:size(values,2)
   values_bits = [values_bits bit_mapping(getindex(values(count),M), :)]; 
end 
values_bits = values_bits';

[x,~] = create_filters(values, tx_filter_type, "RC", Nsamples, Ns, Nsymbols, Fsim, Rs, Rs, Freq, roll_off);

%% Eb/No penalty
    
    f3dB_coef = [0.4 0.5 0.8 1 1.5 2];
    f3dB_vect = f3dB_coef * Rs;
    penalty = zeros(size(f3dB_vect, 2),1);
    figure;

    
    dim = 1;
    
    % The simulation is run for different values of bandwidth
    for f3dB = f3dB_vect

    for ii = 1:length(EbNo_dB)
        %% AWGN channel 
        Ps = var(x);            % Signal Power
        No = (Ps/Rb)./EbNo(ii);
        Pn = No/2*Fsim;
        noise = sqrt(Pn).*randn(Nsamples,1); % noise signal
        y = x+noise; % add WGN to signal
        
        
        %% Apply filter
        % In frequency domain
        [~, H] = create_filters(values, tx_filter_type, "RC", Nsamples, Ns, Nsymbols, Fsim, Rs, f3dB, Freq, roll_off);
        Y = fftshift(fft(y));
        R = H.*Y;
        r = ifft(ifftshift(R));

        %% Ber counting and optimum sampling time
        
        % Get symbols from signal amplitude
        min_err = 1e10;
        
        for k = 1:Ns
            rk = r(k:Ns:end);
            rkd = decode_values(rk,M);
            err(k) = sum(rkd~=values);
            if err(k) < min_err
               min_err = err(k);
               rkd_opt = rkd;
            end
        end

        % Find BER
        errors = 0;
        for count = 1:size(rkd_opt,2)
           if rkd_opt(count) ~= values(count)
               % if the symbol is wrong, the Hamming distance between the two
               % symbols is computed
               errors = errors + sum(bit_mapping(getindex(rkd_opt(count),M),:) ~= bit_mapping(getindex(values(count),M),:));
           end
        end

        BER(ii)  = errors / Nbits;

    end 
        semilogy(EbNo_dB, BER)
        grid on
        hold on
        % Define the value without the penalty of the RC filter
        % (values are found experimentally from BER plots)
        if M == 4
            no_penalty_value = 11;
        else 
            no_penalty_value = 15;
        end
        
        try 
            penalty(dim) = abs(interp1(log10(BER(BER ~= 0 & BER ~= 0.0001)), EbNo_dB(BER ~= 0 & BER ~= 0.0001), -3) - no_penalty_value);
            dim = dim + 1;
        catch
            penalty(dim) = NaN;
            dim = dim + 1;
        end
        
    end 
    
    legend('0.4', '0.5', '0.8', '1', '1.5', '2');
    xlabel('Eb/No [dB]');
    ylabel('BER');

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

