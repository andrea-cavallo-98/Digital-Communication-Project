%
% Andrea Cavallo
% matricola 245715
%
% PROJECT #2
% Main file
%
% DESCRIPTION:
% This file contains the final plots for each step of the project.
% The data plotted are obtained from functions described in the
% other files, run every time with different parameters.
%


clear all
close all
clc

%% General parameters

Nbits = 3e5;        % Number of bits 
%Nsymbols = 1e5;     % Number of symbols
Rb = 1e9;           % Bit rate
Ns = 8;             % Number of samples
roll_off = 0.5;     % Roll-off for SRRC
f3dB_coeff = 1;     % Coefficient of Rs for RC filter bandwidth

%% 2PAM simulations
% Run this section for the plots in step 1, 2, 3

% NRZ pulse shape
[Freq_nrz, PSDx_nrz, EbNo_dB_nrz, BERth_nrz, BER_nrz] = project_02_2PAM( "NRZ", "NRZ", Nbits, Rb, Ns, roll_off, f3dB_coeff);
% RZ pulse shape
[Freq_rz, PSDx_rz, EbNo_dB_rz, BERth_rz, BER_rz] = project_02_2PAM( "RZ", "RZ", Nbits, Rb, Ns, roll_off, f3dB_coeff);
% NRZ pulse shape RC filter
[Freq_rc_nrz, PSDx_rc_nrz, EbNo_dB_rc_nrz, BERth_rc_nrz, BER_rc_nrz] = project_02_2PAM( "NRZ", "RC", Nbits, Rb, Ns, roll_off, f3dB_coeff );
% RZ pulse shape RC filter
[Freq_rc_rz, PSDx_rc_rz, EbNo_dB_rc_rz, BERth_rc_rz, BER_rc_rz] = project_02_2PAM( "RZ", "RC", Nbits, Rb, Ns, roll_off, f3dB_coeff );
% SRRC
[Freq_srrc_005, PSDx_srrc_005, EbNo_dB_srrc_005, BERth_srrc_005, BER_srrc_005] = project_02_2PAM( "SRRC", "SRRC", Nbits, Rb, Ns, 0.05, f3dB_coeff );
[Freq_srrc_05, PSDx_srrc_05, EbNo_dB_srrc_05, BERth_srrc_05, BER_srrc_05] = project_02_2PAM( "SRRC", "SRRC", Nbits, Rb, Ns, 0.5, f3dB_coeff );
[Freq_srrc_1, PSDx_srrc_1, EbNo_dB_srrc_1, BERth_srrc_1, BER_srrc_1] = project_02_2PAM( "SRRC", "SRRC", Nbits, Rb, Ns, 1, f3dB_coeff );

%% 4PAM simulations
% Run this section for the plots in step 1, 2, 3, 5

% NRZ pulse shape
[Freq_4_nrz, PSDx_4_nrz, EbNo_dB_4_nrz, theoryBer_4_nrz, theorySer_4_nrz, BER_4_nrz, SER_4_nrz] = project_02_MPAM(4, 1, "NRZ", "NRZ", roll_off, Nbits / 2, Ns, Rb, f3dB_coeff);
% RZ pulse shape
[Freq_4_rz, PSDx_4_rz, EbNo_dB_4_rz, theoryBer_4_rz, theorySer_4_rz, BER_4_rz, SER_4_rz] = project_02_MPAM(4, 1, "RZ", "RZ", roll_off, Nbits / 2, Ns, Rb, f3dB_coeff);
% NRZ pulse shape RC filter
[Freq_4_nrz_rc, PSDx_4_nrz_rc, EbNo_dB_4_nrz_rc, theoryBer_4_nrz_rc, theorySer_4_nrz_rc, BER_4_nrz_rc, SER_4_nrz_rc] = project_02_MPAM(4, 1, "NRZ", "RC", roll_off, Nbits / 2, Ns, Rb, f3dB_coeff);
% RZ pulse shape RC filter
[Freq_4_rz_rc, PSDx_4_rz_rc, EbNo_dB_4_rz_rc, theoryBer_4_rz_rc, theorySer_4_rz_rc, BER_4_rz_rc, SER_4_rz_rc] = project_02_MPAM(4, 1, "RZ", "RC", roll_off, Nbits / 2, Ns, Rb, f3dB_coeff);
% SRRC
[Freq_4_srrc_005, PSDx_4_srrc_005, EbNo_dB_4_srrc_005, theoryBer_4_srrc_005, theorySer_4_srrc_005, BER_4_srrc_005, SER_4_srrc_005] = project_02_MPAM(4, 1, "SRRC", "SRRC", 0.05, Nbits / 2, Ns, Rb, f3dB_coeff);
[Freq_4_srrc_05, PSDx_4_srrc_05, EbNo_dB_4_srrc_05, theoryBer_4_srrc_05, theorySer_4_srrc_05, BER_4_srrc_05, SER_4_srrc_05] = project_02_MPAM(4, 1, "SRRC", "SRRC", 0.5, Nbits / 2, Ns, Rb, f3dB_coeff);
[Freq_4_srrc_1, PSDx_4_srrc_1, EbNo_dB_4_srrc_1, theoryBer_4_srrc_1, theorySer_4_srrc_1, BER_4_srrc_1, SER_4_srrc_1] = project_02_MPAM(4, 1, "SRRC", "SRRC", 1, Nbits / 2, Ns, Rb, f3dB_coeff);

%% 8PAM simulations
% Run this section for the plots in step 1, 2, 3, 5

% NRZ pulse shape
[Freq_8_nrz, PSDx_8_nrz, EbNo_dB_8_nrz, theoryBer_8_nrz, theorySer_8_nrz, BER_8_nrz, SER_8_nrz] = project_02_MPAM(8, 1, "NRZ", "NRZ", roll_off, Nbits / 3, Ns, Rb, f3dB_coeff);
% RZ pulse shape
[Freq_8_rz, PSDx_8_rz, EbNo_dB_8_rz, theoryBer_8_rz, theorySer_8_rz, BER_8_rz, SER_8_rz] = project_02_MPAM(8, 1, "RZ", "RZ", roll_off, Nbits / 3, Ns, Rb, f3dB_coeff);
% NRZ pulse shape RC filter
[Freq_8_nrz_rc, PSDx_8_nrz_rc, EbNo_dB_8_nrz_rc, theoryBer_8_nrz_rc, theorySer_8_nrz_rc, BER_8_nrz_rc, SER_8_nrz_rc] = project_02_MPAM(8, 1, "NRZ", "RC", roll_off, Nbits / 3, Ns, Rb, f3dB_coeff);
% RZ pulse shape RC filter
[Freq_8_rz_rc, PSDx_8_rz_rc, EbNo_dB_8_rz_rc, theoryBer_8_rz_rc, theorySer_8_rz_rc, BER_8_rz_rc, SER_8_rz_rc] = project_02_MPAM(8, 1, "RZ", "RC", roll_off, Nbits / 3, Ns, Rb, f3dB_coeff);
% SRRC
[Freq_8_srrc_005, PSDx_8_srrc_005, EbNo_dB_8_srrc_005, theoryBer_8_srrc_005, theorySer_8_srrc_005, BER_8_srrc_005, SER_8_srrc_005] = project_02_MPAM(8, 1, "SRRC", "SRRC", 0.005, Nbits / 3, Ns, Rb, f3dB_coeff);
[Freq_8_srrc_05, PSDx_8_srrc_05, EbNo_dB_8_srrc_05, theoryBer_8_srrc_05, theorySer_8_srrc_05, BER_8_srrc_05, SER_8_srrc_05] = project_02_MPAM(8, 1, "SRRC", "SRRC", 0.5, Nbits / 3, Ns, Rb, f3dB_coeff);
[Freq_8_srrc_1, PSDx_8_srrc_1, EbNo_dB_8_srrc_1, theoryBer_8_srrc_1, theorySer_8_srrc_1, BER_8_srrc_1, SER_8_srrc_1] = project_02_MPAM(8, 1, "SRRC", "SRRC", 1, Nbits / 3, Ns, Rb, f3dB_coeff);


%% Step 1: Plot and compare transmitted spectra

% 2PAM
figure;

plot(Freq_nrz,10*log10(PSDx_nrz),'r')
hold on
plot(Freq_rz,10*log10(PSDx_rz),'b');
hold on
plot(Freq_srrc_05,10*log10(PSDx_srrc_05),'k');
grid on
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC')
title('Spectral Analysis of transmitted signal for 2PAM');

% 4PAM

figure;
plot(Freq_4_nrz,10*log10(PSDx_4_nrz),'r')
hold on
plot(Freq_4_rz,10*log10(PSDx_4_rz),'b');
hold on
plot(Freq_4_srrc_05,10*log10(PSDx_4_srrc_05),'k');
grid on
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC')
title('Spectral Analysis of transmitted signal for 4PAM');

% 8PAM
figure;
plot(Freq_8_nrz,10*log10(PSDx_8_nrz),'r')
hold on
plot(Freq_8_rz,10*log10(PSDx_8_rz),'b');
hold on
plot(Freq_8_srrc_05,10*log10(PSDx_8_srrc_05),'k');
grid on
xlabel('Freq [Hz]');
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC')
title('Spectral Analysis of transmitted signal for 8PAM');

% NRZ filters
figure;
plot(Freq_nrz, 10*log10(PSDx_nrz), 'r');
hold on;
plot(Freq_4_nrz, 10*log10(PSDx_4_nrz), 'b');
hold on;
plot(Freq_8_nrz, 10*log10(PSDx_8_nrz), 'k');
grid on;
xlabel('Freq [Hz]')
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('2PAM', '4PAM', '8PAM');
title('NRZ pulse shape');
axis([-3e9 3e9 -45 0])

% RZ filters
figure;
plot(Freq_rz, 10*log10(PSDx_rz), 'r');
hold on;
plot(Freq_4_rz, 10*log10(PSDx_4_rz), 'b');
hold on;
plot(Freq_8_rz, 10*log10(PSDx_8_rz), 'k');
grid on;
xlabel('Freq [Hz]')
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('2PAM', '4PAM', '8PAM');
title('RZ pulse shape');
axis([-3e9 3e9 -35 0])

% SRRC filters
figure;
plot(Freq_srrc_005, 10*log10(PSDx_srrc_005), 'r');
hold on;
plot(Freq_4_srrc_005, 10*log10(PSDx_4_srrc_005), 'b');
hold on;
plot(Freq_8_srrc_005, 10*log10(PSDx_8_srrc_005), 'k');
grid on;
xlabel('Freq [Hz]')
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('2PAM', '4PAM', '8PAM');
title('SRRC pulse shape with roll-off = 0.05');
axis([-3e9 3e9 -45 0])

figure;
plot(Freq_srrc_1, 10*log10(PSDx_srrc_1), 'r');
hold on;
plot(Freq_4_srrc_1, 10*log10(PSDx_4_srrc_1), 'b');
hold on;
plot(Freq_8_srrc_1, 10*log10(PSDx_8_srrc_1), 'k');
grid on;
xlabel('Freq [Hz]')
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('2PAM', '4PAM', '8PAM');
title('SRRC pulse shape with roll-off = 1');
axis([-3e9 3e9 -45 0])

figure;
plot(Freq_srrc_05, 10*log10(PSDx_srrc_05), 'r');
hold on;
plot(Freq_4_srrc_05, 10*log10(PSDx_4_srrc_05), 'b');
hold on;
plot(Freq_8_srrc_05, 10*log10(PSDx_8_srrc_05), 'k');
grid on;
xlabel('Freq [Hz]')
ylabel('Power Spectral Density [dB/sqrt(Hz)]')
legend('2PAM', '4PAM', '8PAM');
title('SRRC pulse shape with roll-off = 0.5');
axis([-3e9 3e9 -45 0])

%% Step 2: Transmission in AWGN channel with matched filter

% 2PAM
figure;
semilogy(EbNo_dB_nrz, BER_nrz, 'Marker','o','MarkerEdgeColor','r','MarkerSize',7,'LineStyle','none');
hold on
semilogy(EbNo_dB_rz, BER_rz, 'Marker','x','MarkerEdgeColor','b','MarkerSize',10,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_srrc_05, BER_srrc_05, 'Marker','d','MarkerEdgeColor','g','MarkerSize',5,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_nrz, BERth_nrz, 'k-');
grid on
xlabel('EbNo [dB]')
ylabel('BER')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC', 'Theoretical curve');
title('2PAM results');

% 4PAM
figure;
semilogy(EbNo_dB_4_nrz, BER_4_nrz, 'Marker','o','MarkerEdgeColor','r','MarkerSize',7,'LineStyle','none');
hold on
semilogy(EbNo_dB_4_rz, BER_4_rz, 'Marker','x','MarkerEdgeColor','b','MarkerSize',10,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_4_srrc_05, BER_4_srrc_05, 'Marker','d','MarkerEdgeColor','g','MarkerSize',5,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_4_nrz, theoryBer_4_nrz, 'k-');
grid on
xlabel('EbNo [dB]')
ylabel('BER')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC', 'Theoretical curve');
title('4PAM BER results');

figure;
semilogy(EbNo_dB_4_nrz, SER_4_nrz, 'Marker','o','MarkerEdgeColor','r','MarkerSize',7,'LineStyle','none');
hold on
semilogy(EbNo_dB_4_rz, SER_4_rz, 'Marker','x','MarkerEdgeColor','b','MarkerSize',10,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_4_srrc_05, SER_4_srrc_05, 'Marker','d','MarkerEdgeColor','g','MarkerSize',5,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_4_nrz, theorySer_4_nrz, 'k-');
grid on
xlabel('EbNo [dB]')
ylabel('SER')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC', 'Theoretical curve');
title('4PAM SER results');


% 8PAM
figure;
semilogy(EbNo_dB_8_nrz, BER_8_nrz, 'Marker','o','MarkerEdgeColor','r','MarkerSize',7,'LineStyle','none');
hold on
semilogy(EbNo_dB_8_rz, BER_8_rz, 'Marker','x','MarkerEdgeColor','b','MarkerSize',10,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_8_srrc_05, BER_8_srrc_05, 'Marker','d','MarkerEdgeColor','g','MarkerSize',5,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_8_nrz, theoryBer_8_nrz, 'k-');
grid on
xlabel('EbNo [dB]')
ylabel('BER')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC', 'Theoretical curve');
title('8PAM BER results');

figure;
semilogy(EbNo_dB_8_nrz, SER_8_nrz, 'Marker','o','MarkerEdgeColor','r','MarkerSize',7,'LineStyle','none');
hold on
semilogy(EbNo_dB_8_rz, SER_8_rz, 'Marker','x','MarkerEdgeColor','b','MarkerSize',10,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_8_srrc_05, SER_8_srrc_05, 'Marker','d','MarkerEdgeColor','g','MarkerSize',5,'LineStyle','none', 'LineWidth', 1);
hold on
semilogy(EbNo_dB_8_nrz, theorySer_8_nrz, 'k-');
grid on
xlabel('EbNo [dB]')
ylabel('SER')
legend('NRZ pulse shape', 'RZ pulse shape', 'SRRC', 'Theoretical curve');
title('8PAM SER results');

% Different channels comparison BER
figure;
semilogy(EbNo_dB_8_nrz, BER_8_nrz, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_4_nrz, BER_4_nrz, 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_nrz, BER_nrz, 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_8_nrz, theoryBer_8_nrz, 'r-');
hold on
semilogy(EbNo_dB_4_nrz, theoryBer_4_nrz, 'b-');
hold on
semilogy(EbNo_dB_nrz, BERth_nrz, 'k-');
hold on
grid on
xlabel('Eb/No [dB]');
ylabel('BER')
legend('8PAM','4PAM','2PAM','8PAM theory', '4PAM theory', '2PAM theory');
title('BER for NRZ pulse shape with matched filter');

% Different channels comparison SER
figure;
semilogy(EbNo_dB_8_nrz, SER_8_nrz, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_4_nrz, SER_4_nrz, 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_nrz, BER_nrz, 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_8_nrz, theorySer_8_nrz, 'r-');
hold on
semilogy(EbNo_dB_4_nrz, theorySer_4_nrz, 'b-');
hold on
semilogy(EbNo_dB_nrz, BERth_nrz, 'k-');
hold on
grid on
xlabel('Eb/No [dB]');
ylabel('SER')
legend('8PAM','4PAM','2PAM','8PAM theory', '4PAM theory', '2PAM theory');
title('SER for NRZ pulse shape with matched filter');

%% Step 3: Simulations

% NOTE: The simulation may stop in error because of the function
% interp1 used in project_02_2PAM_RCpenalty and project_02_MPAM_RCpenalty.
% If it gives an error, it is sufficient to run it again. 

% 2PAM
[penalty_2_NRZ, f3dB_vect_2_NRZ] = project_02_2PAM_RCpenalty("NRZ", roll_off, Nbits, Ns, Rb);
title('2PAM BER for different RC filters with NRZ pulse shape');
[penalty_2_RZ, f3dB_vect_2_RZ] = project_02_2PAM_RCpenalty("RZ", roll_off, Nbits, Ns, Rb);
title('2PAM BER for different RC filters with RZ pulse shape');

% 4PAM
[penalty_4_NRZ, f3dB_vect_4_NRZ] = project_02_MPAM_RCpenalty(4, 1, "NRZ", roll_off, Nbits / 2, Ns, Rb);
title('4PAM BER for different RC filters with NRZ pulse shape');
[penalty_4_RZ, f3dB_vect_4_RZ] = project_02_MPAM_RCpenalty(4, 1, "RZ", roll_off, Nbits / 2, Ns, Rb);
title('4PAM BER for different RC filters with RZ pulse shape');

% 8PAM
[penalty_8_NRZ, f3dB_vect_8_NRZ] = project_02_MPAM_RCpenalty(8, 1, "NRZ", roll_off, Nbits / 3, Ns, Rb);
title('8PAM BER for different RC filters with NRZ pulse shape');
[penalty_8_RZ, f3dB_vect_8_RZ] = project_02_MPAM_RCpenalty(8, 1, "RZ", roll_off, Nbits / 3, Ns, Rb);
title('8PAM BER for different RC filters with RZ pulse shape');


%% Step 3: Transmission in AWGN channel with RC filter

% Penalty with NRZ
figure;
plot(f3dB_vect_2_NRZ, penalty_2_NRZ,'r', 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
hold on
plot(f3dB_vect_4_NRZ, penalty_4_NRZ,'b', 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
hold on
plot(f3dB_vect_8_NRZ, penalty_8_NRZ,'k', 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
grid on
xlabel('3dB bandwidth of RC filter / Rs');
ylabel('Eb/No penalty [dB]');
title('Eb/No penalty for NRZ pulse shape');
legend('2PAM', '4PAM', '8PAM');


% Penalty with RZ
figure;
plot(f3dB_vect_2_RZ, penalty_2_RZ, 'r', 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
hold on
plot(f3dB_vect_4_RZ, penalty_4_RZ, 'b', 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
hold on
plot(f3dB_vect_8_RZ, penalty_8_RZ, 'k', 'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
grid on
xlabel('3dB bandwidth of RC filter / Rs');
ylabel('Eb/No penalty [dB]');
title('Eb/No penalty for RZ pulse shape');
legend('2PAM', '4PAM', '8PAM');

%%% Spectral analysis

% 2PAM
figure;
H2_rz = 1./(1+1i*Freq_nrz./(0.7*Rb));
H2_nrz = 1./(1+1i*Freq_nrz./(0.5*Rb));
H2_rz = H2_rz./max(abs(H2_rz));
H2_nrz = H2_nrz./max(abs(H2_nrz));
plot(Freq_nrz, (PSDx_nrz), 'b', 'LineWidth', 1);
hold on
plot(Freq_nrz, (PSDx_rz), 'r', 'LineWidth', 1);
hold on
plot(Freq_nrz, (abs(H2_rz)), 'r--', 'LineWidth', 1);
hold on
plot(Freq_nrz, (abs(H2_nrz)), 'b--', 'LineWidth', 1);
grid on
xlabel('Freq [Hz]');
ylabel('Power spectral density')
legend('NRZ pulse shape', 'RZ pulse shape', 'Best RC filter (RZ)', 'Best RC filter (NRZ)');
title('Spectral analysis for 2PAM');

% 4PAM
figure;
H4_rz = 1./(1+1i*Freq_4_nrz./(1*Rb/2));
H4_nrz = 1./(1+1i*Freq_4_nrz./(0.5*Rb/2));
H4_rz = H4_rz./max(abs(H4_rz));
H4_nrz = H4_nrz./max(abs(H4_nrz));
plot(Freq_4_nrz, (PSDx_4_nrz), 'b', 'LineWidth', 1);
hold on
plot(Freq_4_nrz, (PSDx_4_rz), 'r', 'LineWidth', 1);
hold on
plot(Freq_4_nrz, (abs(H4_nrz)), 'b--', 'LineWidth', 1);
hold on
plot(Freq_4_nrz, (abs(H4_rz)), 'r--', 'LineWidth', 1);
grid on
xlabel('Freq [Hz]');
ylabel('Power spectral density')
legend('NRZ pulse shape', 'RZ pulse shape', 'Best RC filter (NRZ)', 'Best RC filter (RZ)');
title('Spectral analysis for 4PAM');

% 8PAM
figure;
H8_rz = 1./(1+1i*Freq_8_nrz./(1*Rb/3));
H8_nrz = 1./(1+1i*Freq_8_nrz./(0.8*Rb/3));
H8_rz = H8_rz./max(abs(H8_rz));
H8_nrz = H8_nrz./max(abs(H8_nrz));
plot(Freq_8_nrz, (PSDx_8_nrz), 'b', 'LineWidth', 1);
hold on
plot(Freq_8_nrz, (PSDx_8_rz), 'r', 'LineWidth', 1);
hold on
plot(Freq_8_nrz, abs(H8_nrz), 'b--', 'LineWidth', 1);
hold on
plot(Freq_8_nrz, abs(H8_rz), 'r--', 'LineWidth', 1);
grid on
xlabel('Freq [Hz]');
ylabel('Power spectral density')
legend('NRZ pulse shape', 'RZ pulse shape', 'Best RC filter (NRZ)', 'Best RC filter (RZ)');
title('Spectral analysis for 8PAM');

%% BER with ISI: simulation
% Run this section for the plots in step 4
[EbNo_dB_ISI_nrz, BER_mean_ISI_nrz, BERth_ISI_nrz, BER_worst_ISI_nrz, BER_ISI_nrz] = project_02_2PAM_ISI( "NRZ", "RC", Nbits, Rb, 10, roll_off, 0.75 );
title('Eye diagram with NRZ pulse shape');
xlabel('Time');
ylabel('Amplitude');
[EbNo_dB_ISI_rz, BER_mean_ISI_rz, BERth_ISI_rz, BER_worst_ISI_rz, BER_ISI_rz] = project_02_2PAM_ISI( "RZ", "RC", Nbits, Rb, 10, roll_off, 0.75 );
title('Eye diagram with RZ pulse shape');
xlabel('Time');
ylabel('Amplitude');

%% Step 4: BER with ISI by means of the eye diagram

% NRZ
figure;
semilogy(EbNo_dB_ISI_nrz, BER_ISI_nrz, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_ISI_nrz, BER_worst_ISI_nrz, 'b-');
hold on
semilogy(EbNo_dB_ISI_nrz, BER_mean_ISI_nrz, 'k-');
hold on
semilogy(EbNo_dB_ISI_nrz, BERth_ISI_nrz, 'k--');
xlabel('Eb/No [dB]');
ylabel('BER');
legend('Experimental values','Theoretical value using worst cases', 'Theoretical value using mean values', 'Theoretical value without ISI');
grid on
title('BER results with ISI for NRZ')

% NRZ
figure;
semilogy(EbNo_dB_ISI_rz, BER_ISI_rz, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_ISI_rz, BER_worst_ISI_rz, 'b-');
hold on
semilogy(EbNo_dB_ISI_rz, BER_mean_ISI_rz, 'k-');
hold on
semilogy(EbNo_dB_ISI_rz, BERth_ISI_rz, 'k--');
xlabel('Eb/No [dB]');
ylabel('BER');
legend('Experimental values','Theoretical value using worst cases', 'Theoretical value using mean values', 'Theoretical values without ISI');
grid on
title('BER results with ISI for RZ')

%% Non-Gray coding: simulation
% Run this section for the plots in step 5

% 4PAM NRZ pulse shape
[Freq_ng_4, PSDx_ng_4, EbNo_dB_ng_4, theoryBer_ng_4, theorySer_ng_4, BER_ng_4, SER_ng_4] = project_02_MPAM(4, 0, "NRZ", "NRZ", roll_off, Nbits / 2, Ns, Rb, f3dB_coeff);

% 8PAM NRZ pulse shape
[Freq_ng_8, PSDx_ng_8, EbNo_dB_ng_8, theoryBer_ng_8, theorySer_ng_8, BER_ng_8, SER_ng_8] = project_02_MPAM(8, 0, "NRZ", "NRZ", roll_off, Nbits / 3, Ns, Rb, f3dB_coeff);

%% Step 5: Non-Gray coding: plots

% 4PAM
figure;
semilogy(EbNo_dB_ng_4, BER_ng_4, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_4_nrz, BER_4_nrz, 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_ng_4, theoryBer_ng_4, 'r--');
hold on
semilogy(EbNo_dB_4_nrz, theoryBer_4_nrz, 'b--');
grid on
xlabel('EbNo [dB]')
ylabel('BER')
legend('NRZ pulse shape (non-Gray coding)', 'NRZ pulse shape (Gray coding)', 'Theoretical BER (non-Gray coding)', 'Theoretical BER (Gray coding)');
title('4PAM BER with non-Gray coding');

% 8PAM
figure;
semilogy(EbNo_dB_ng_8, BER_ng_8, 'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_8_nrz, BER_8_nrz, 'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5,'LineStyle','none');
hold on
semilogy(EbNo_dB_ng_8, theoryBer_ng_8, 'r--');
hold on
semilogy(EbNo_dB_8_nrz, theoryBer_8_nrz, 'b--');
grid on
xlabel('EbNo [dB]')
ylabel('BER')
legend('NRZ pulse shape (non-Gray coding)', 'NRZ pulse shape (Gray coding)','Theoretical BER (non-Gray coding)', 'Theoretical BER (Gray coding)');
title('8PAM BER with non-Gray coding');
