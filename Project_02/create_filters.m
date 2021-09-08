%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #2
% Create filters
%
% DESCRIPTION:
% This file contains the functions to generate and impress the pulse shape
% to the transmitter signal and also to generate the transfer function for
% the receiver filter
%


function [x, H] = create_filters(values, tx_filter_type, rx_filter_type, Nsamples, Ns, Nbits, Fsim, Rs, f3dB, Freq, roll_off)

% Transmitter filter
if tx_filter_type == "NRZ"
    s = ones(Ns,1);
    x = zeros(Nsamples,1);
    for n = 1:Nbits
       idx = (n-1)*Ns +(1:Ns);
      x(idx,1) = values(n).*s;
    end
    
end
if tx_filter_type == "RZ"
    s = [ones(Ns/2,1); zeros(Ns/2,1)];
    x = zeros(Nsamples,1);
    for n = 1:Nbits
        idx = (n-1)*Ns +(1:Ns);
        x(idx,1) = values(n).*s;
    end
end
if tx_filter_type == "SRRC"
    Symbols = values;
    f = [-Nsamples/2:Nsamples/2-1].*Fsim/Nsamples;
    TF = zeros(Nsamples, 1);
    alpha = 0.5; 
    a = roll_off;
    W = Rs / 2;
    f2 = (1+a)*W;
    f1 = (1-a)*W;
    idx1=find(abs(f)<=f1);
    idx2=find(abs(f)>=f1 & abs(f) <= f2);
    TF(idx2) = 0.5 * (1-sin(pi*(abs(f(idx2))-W)/(2*a*W)));
    TF(idx1) = 1;
    TF = TF.^alpha;
    TF = TF.*Ns;

    x_tx = upsample(Symbols,Ns); % train of deltas
    x = ifft(fftshift(TF)'.*fft(x_tx)); % convolution between train of deltas and filter
    x = x';
end


% Receiver filter
if rx_filter_type == "NRZ"
    H = (1./Ns).*conj(fftshift(fft(s,Nsamples)));
end
if rx_filter_type == "RZ"
    H = (2./Ns).*conj(fftshift(fft(s,Nsamples)));
end
if rx_filter_type == "RC"
    H = 1./(1+1i*Freq./f3dB);
    H = H./max(abs(H));
end
if rx_filter_type == "SRRC"
    H = TF/Ns;
end

end