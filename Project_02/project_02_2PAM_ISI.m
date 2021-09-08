%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #2
% 2-PAM system with ISI
%
% DESCRIPTION:
% This file performs the evaluation of ISI by means of eye diagram
% analysis. Plots are generated for the eye diagrams.
%



function [EbNo_dB, BER_mean, BERth, BER_worst, BER] = project_02_2PAM_ISI( tx_filter_type, rx_filter_type, Nbits, Rb, Ns, roll_off, f3dB_coeff )


%% General parameters 

Nsamples = Ns*Nbits;                % Number of samples
Fsim = Ns*Rb;                       % Simulation Bandwidth [Hz]
EbNo_dB = [2:12];                    % Eb/No [dB]
EbNo = 10.^(EbNo_dB*0.1);           % Eb/No
Rs = Rb;                            % Symbol rate
f3dB = f3dB_coeff * Rs;             % Bandwidth of RC filter
nfft = Nsamples;                    % Samples for Fourier transform
stepFreq = Fsim/nfft;               % Step-frequency
maxFreq = +Fsim/2-stepFreq;         % Max frequency
minFreq = -Fsim/2;                  % Min frequency
Freq = (minFreq:stepFreq:maxFreq)'; % Frequency vector

%% Transmitter

% Bit generation 
Bits = randi([0 1], Nbits, 1);
% Antipodal representation
values(Bits == 0) = -1;
values(Bits == 1) = 1;

[x,H] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nbits, Fsim, Rs, f3dB, Freq, roll_off);

%% Plot eye diagram

% Apply filter
X = fftshift(fft(x));
R = H.*X;
r = real(ifft(ifftshift(R)));

% Plot
lines = eyePlot(r,1000, 2*Ns,floor(Ns/2),2);


%% Find BER with ISI using formulas

% Find the optimal sampling instant
maxdiff = 0;
for count = 1:Ns
    
    % Minimum positive value and maximum negative value of eye
    min_pos = min(lines(lines(:,count) > 0,count));
    max_neg = max(lines(lines(:,count) < 0,count));
    diff = abs(min_pos - max_neg);
    % Identify the sampling instant of maximum eye opening
    if diff > maxdiff
        kopt = count;
        maxdiff = diff;
    end
    
end

% First method: find upper and lower averages
up_mean = mean(lines(lines(:,kopt)>0,kopt));
down_mean = mean(lines(lines(:,kopt)<=0,kopt));

% Second method: worst cases
up_worst = min(lines(lines(:,kopt)>0,kopt));
down_worst = max(lines(lines(:,kopt)<=0,kopt));

% Equivalent bandwidth
Beq = sum(abs(H).^2)/2.*stepFreq;

% Noise power and BER
A = var(x);
k = (up_mean - down_mean)^2/8*Rb./(A*Beq);
k_worst = (up_worst - down_worst)^2/8*Rb./(A*Beq);
BER_mean = 0.5*erfc(sqrt(k.*EbNo));
BER_worst = 0.5*erfc(sqrt(k_worst.*EbNo));
BERth = 0.5*erfc(sqrt(EbNo));

%% Transmission and filtering (simulation)

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


%% Functions

% Plot the eye diagram and return the values
function [lines] = eyePlot(x,Ntraces,SpT,Offset,Period)
time_plot = linspace(-Period/2,Period/2,SpT);
x = x(Offset:end,:);
figure()
lines = zeros(Ntraces,size((1:SpT),2));
for n = 1:Ntraces
    idx = (n-1)*SpT+(1:SpT);
    lines(n,:) = x(idx);
    plot(time_plot,x(idx), 'b-');
    hold on;
end
end
