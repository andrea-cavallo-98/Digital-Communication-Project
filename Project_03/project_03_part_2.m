%
% Andrea Cavallo
% matricola 245715
%
%
% PROJECT #3
% PART 2: RF Spectrum Analyzer using
% RTL-SDR receiver
%



clear all
close all
clc

%% General parameters

LeftFrequency = 174e6; % Left frequency of the sweep
RightFrequency = 210e6; % Right frequency of the sweep
SampleRate = 1.2e6; % Rate of the acquisition (partial overlap)
EffectiveSampleRate = 1.0e6; % Bandwidth in the plot (overlaps removed)
NsBlock = 1024;
Nblocks = 1000;

fm_radios = [88.2 88.5 88.7 89 89.3 89.3 89.7 90 90.3 90.7 90.9 91.2 91.5 91.8 92.1 92.4 ...
    92.7 93 93.3 93.6 93.9 94.2 94.4 94.7 95 95.3 95.6 95.9 96.2 96.4 96.7 97 97.3 97.6 97.9 ...
    98.2 98.5 98.7 99 99.3 99.6 99.9 100.2 100.5 100.75 101 101.25 101.5 101.8 ...
    102.1 102.5 102.8 103.1 103.3 103.5 103.7 104 104.2 104.45 104.6 105 105.25 105.5 105.75 ...
    106 106.3 106.6 106.9 107.4 107.7];

dab_radios = [215.072 223.936 225.648 227.360 229.072];

% dab_radios = [174.928 176.640 178.352 180.064 181.936 183.648 185.360 187.072 188.928...
% 190.640 192.352 194.064 195.936 197.648 199.360 201.072 202.928 204.640...
% 206.352 208.064 209.936 211.648 213.360 215.072 216.928 218.640 220.352...
% 222.064 223.936 225.648 227.360 229.072 230.784 232.496 234.208 235.776...
% 237.488 239.200 ];

%% Simulation

Frequencies = LeftFrequency:EffectiveSampleRate:RightFrequency;

figure;
rep = 0;

for CenterFrequency = Frequencies


    %% Signal acquisition

    % Acquire the signal
    hRadio = comm.SDRRTLReceiver('CenterFrequency', CenterFrequency, ...
        'SampleRate', SampleRate, 'EnableTunerAGC', true, ...
        'SamplesPerFrame', NsBlock, 'OutputDataType', 'single');

    % Extract signal
    if ~isempty(sdrinfo(hRadio.RadioAddress))
        x_saved = NaN*ones(Nblocks*NsBlock, 1);
        for Counter = 1:Nblocks
            [x, len, lost(Counter)] = step(hRadio);
            x_saved((Counter-1)*NsBlock+1:(Counter-1)*NsBlock+NsBlock) = x;
        end
    else
        warning('SDR Device not connected');
    end
    release(hRadio);

    %% Spectral analysis

    % Apply Bartlett periodogram
    [f, X] = project_03_myBartlett2(x_saved, round(size(x_saved, 1) / 100), SampleRate, LeftFrequency + EffectiveSampleRate * rep);
    
    % Eliminate borders to avoid superimposition
    to_cut = abs(SampleRate - EffectiveSampleRate) / 2;
    idx = round(to_cut/SampleRate * round(size(x_saved, 1) / 100)) : round((SampleRate - to_cut)/SampleRate * round(size(x_saved, 1) / 100));
    
    % Plot
    plot(f(idx) / 1e6, 10*log10(abs(X(idx))), 'b', 'LineWidth', 1);
    hold on
    
    rep = rep + 1;
    
end

xlabel('Freq [MHz]')
ylabel('Power Spectral Density [dB]');

title('Spectrum Sweep 88 MHz - 108 MHz, Torino');

%% Add expected FM radio frequencies

for fm = fm_radios
plot([fm fm], [-100 100], '--');
hold on
end

%% Add expected DAB radio frequencies

for dab = dab_radios
plot([dab dab], [-100 100], '--');
hold on
end

