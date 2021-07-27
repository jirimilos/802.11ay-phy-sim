%% Performance Analysis of IEEE 802.11 WLAN Technologies
% main simulation batch file; modifications are possible
%
% Authors:	Jiri Milos, Jiri Blumenstein and Ladislav Polak, DREL FEEC BUT, 2020
%
close all
clear all
% clearvars;
clear global
clc

global DEBUG_MODE;
DEBUG_MODE = false;

saveResults = 1;
if saveResults == false
    disp(' ')
    disp('######## RESULTS ARE NOT STORED NOW ########');
    disp(' '), pause(2);
end
%% Define simulation parameters ----------------------------------------
wifi_standard = '802dot11ay'; % 802.11ay physical layer (different PHY layer modes are defined by MCS value)

ChannelBandwidth = [2640]; % [MHz]

UserMode = 'single'; % should be adjustable SU/MU
antMode = 'SISO'; % SISO, RxD, TxD, MIMO
TxRxAnt = 11;


useNUC = false; % use Non-Uniform constellation
use8PSK = false; % use pi/2-8-PSK modulation for MCS12 and MCS13

% GuardInterval = 'normal'; % adjustable {'short', 'normal', 'long'}; normal for now
GuardInterval = 'normal'; % adjustable {'short', 'normal', 'long'}; normal for now
GuardIntervalType = 'GI'; % 'GI' = Guard Interval (Golay sequence), 'CP' = Cyclic Prefix, 'zeros' = zero vector is used

LDPCmatrix = 'normal'; % 'normal' or 'lifted'=long; normal for now
LDPCmatrix = 'lifted'; 

decision_type = 'LLR'; %'HD' - Hard decision values; 'LLR' - Log-likelihood ratio(LLR);

% LENGTH = 300; % user data octet length (1--262143)
LENGTHatMCS = [250,...
    250,300,250,250,250,250,275,250,275,275,275,200,200,200,200,200,200,250,275,275,275]; % user data octet length (1--262143)

% MCSvec = [1:1]; % individual MCS types or string 'all', 'Ctrl', 'all-SC', 'all-OFDM', 'all-LPSC'
% MCSvec = [12 13]; % individual MCS types or string 'all', 'Ctrl', 'all-SC', 'all-OFDM', 'all-LPSC'
% MCSvec = [6,11,16,21]; % individual MCS types or string 'all', 'Ctrl', 'all-SC', 'all-OFDM', 'all-LPSC'
MCSvec = 17:21; % individual MCS types or string 'all', 'Ctrl', 'all-SC', 'all-OFDM', 'all-LPSC'

channelType = 'awgn'; % 'awgn', 'fad' or 'fad_meas' or 'los'
channelKnowledge = 'perfect'; % perfect or estimation

SNR = 10:0.5:30;

% SNR = [0:8, 8.5:0.5:16, 17:20];

N_frames = 200; % number of frames within inner loop (for each SNR value)
n_users = [1 1]; % [AP STA] tx always set to 1

% necessary MCS definition and MCS vector indices definition
mcs_definition;

for i_mcs = i_mcs_vec_aux % mcs loop
    
%     LENGTH = LENGTHatMCS(i_mcs)
    LENGTH = 1000;
    tic;
    % simulation parameters
    load_wifi_params;
    
    % preallocate result variables
    results = result_allocation(SNR,N_frames,wifi_params);
    % main simulation file
    sim_file;
    
    % save results for each MCS
    results = result_all(results);
    % generate results filename
    output_filename = load_filename(SNR,N_frames,i_mcs,wifi_params,wifi_standard, ChannelBandwidth, wifi_params.MCS(i_mcs).coding_type, LENGTH, GuardInterval, decision_type);
    filename_suffix = [];
    if saveResults == 1
        save(fullfile('.\results\nonNUC_PNbothSides\',[output_filename, filename_suffix, '.mat']));
    end
    t_elapsed(1, i_mcs) = toc;
    
end % mcs loop

disp('-------------------------------------------------')
disp(['Elapsed time: ', num2str(sum(t_elapsed)),' seconds.'])
disp(' ')
disp('SNR: ')
disp(SNR)
disp('BER CODED DATA: ')
disp(results.ber_coded_data)
disp('BER UN-CODED DATA: ')
disp(results.ber_uncoded_data)

disp('THR CODED DATA: ')
disp(results.throughput_coded_data)

%% present BER results for single MCS value
if length(MCSvec) == 1
    figure(1)
    subplot(131)
    semilogy(SNR, results.ber_PSDU)
    hold on, grid on
%     semilogy(SNR, results.ber_PSDU_padded)
%     semilogy(SNR, results.ber_coded_data,'k:')
    semilogy(SNR, results.ber_uncoded_data,'k--')
    xlabel('SNR [dB]'), xlim([SNR(1) SNR(end)])
    ylabel('BER [-]'), ylim([1e-6 1])
%     legend('PSDU','PSDU + padding','coded','uncoded','Location','southeast')
    legend('PSDU','PSDU + padding','uncoded','Location','southeast')
    
    subplot(132)
    semilogy(SNR, results.fer_frames_PSDU)
    hold on, grid on
    semilogy(SNR, results.fer_frames_coded_data)
    xlabel('SNR [dB]'), xlim([SNR(1) SNR(end)])
    ylabel('FER [-]'), ylim([1e-6 1])
    legend('PSDU','PSDU + padding')
    
    subplot(133)
    plot(SNR, results.throughput_coded_data,'s-')
    hold on, grid on
    xlabel('SNR [dB]'), xlim([SNR(1) SNR(end)])
    ylabel('Throughput [b/s]')
    xlim([SNR(1), SNR(end)])
%     ylim([0, (floor(ad_MCS_data_rates_vec(i_mcs)/1e9)+1)*1e9])
end