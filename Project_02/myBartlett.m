function [f_Bartlett, F_Bartlett] = myBartlett(Signal, Nfft, Bsim)

if size(Signal, 2) == 1
    Signal = Signal';
end

Repetition = floor(length(Signal)/Nfft);

F_Bartlett = zeros(1, Nfft);
for Counter = 1 : Repetition
    F = fft(Signal((Counter - 1)*Nfft+1 : (Counter - 1)*Nfft + Nfft));
    F = 1/Nfft^2 * abs(fftshift(F)).^2;
    F_Bartlett = F_Bartlett + F;
end
F_Bartlett = F_Bartlett / Repetition;

Df = Bsim/Nfft;
f_Bartlett = [-Bsim/2:Df:Bsim/2-Df];

F_Bartlett = F_Bartlett / max(F_Bartlett);

return