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


function [penalty, f3dB_coef] = project_02_2PAM_RCpenalty(tx_filter_type, roll_off, Nsymbols, Ns, Rb)


%% General parameters 

rx_filter_type="RC";
Nbits=1e5;
Nsamples = Ns*Nbits;             % Number of samples
Fsim = Ns*Rb;                    % Simulation Bandwidth [Hz]
EbNo_dB = [2:15];                % Eb/No [dB]
EbNo = 10.^(EbNo_dB*0.1);        % Eb/No
Rs = Rb;                         % Symbol Rate
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

[x,~] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nbits, Fsim, Rs, Rs, Freq, roll_off);



%% Eb/No penalty

    f3dB_coef = [0.2 0.5 0.7 1 1.3 1.6 2];
    f3dB_vect = f3dB_coef * Rs;
    penalty = zeros(size(f3dB_vect, 2),1);

    dim = 1;
    figure;

    % The simulation of the system is run for different values of
    % the filter bandwidth
    for f3dB = f3dB_vect 

    for ii = 1:length(EbNo_dB)
        
        %% AWGN channel 
        Ps = var(x);            % Signal Power
        No = (Ps/Rb)./EbNo(ii);
        Pn = No/2*Fsim;
        noise = sqrt(Pn).*randn(Nsamples,1); % noise signal
        y = x+noise; % add WGN to the signal
        
    %% Apply filter
        
        % in frequency domain
        [~, H] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nbits, Fsim, Rs, f3dB, Freq);
        Y = fftshift(fft(y));
        R = H.*Y;
        r = ifft(ifftshift(R));

        %% Ber counting and optimum sampling time
        
        Vth = 0.0;
        
        for k = 1:Ns
            rk = r(k:Ns:end);
            rkd = (rk>Vth);
            err(k) =sum(rkd~=Bits);
        end

        [BER(ii),kopt]= min(err./Nbits);


    end 
        semilogy(EbNo_dB, BER)
        grid on
        hold on

        penalty(dim) = abs(interp1(log10(BER(BER ~= 0 & BER ~= 0.0001)), EbNo_dB(BER ~= 0 & BER ~= 0.0001), -3) - 10*log10(erfcinv(1e-3*2)^2));
        dim = dim + 1;

    end 
    legend('0.2', '0.5', '0.7', '1', '1.3', '1.6', '2');
    xlabel('Eb/No [dB]');
    ylabel('BER');
end

