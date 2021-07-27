function [rxObj] = WIFI_RX_ay(wifi_params, rxObj, chanObj, i_mcs, SNR_dB, txObj)
% function rxObj = WIFI_RX_ay(wifi_params,rxObj,channelObj, ChanMod_output, mcs_num, SNR_dB, txObj, noiseVar, var_s)

% WIFI transmitter function according to IEEE 802.11ad standard
% Author:   Jiri Milos, DREL FEEC BUT, 2020

% rxObj = WIFI_RX_ay(wifi_params,rxObj,chanOut, chanMod, i_mcs, SNR(i_snr), txObj, var_txOutput_x_s);

%% DEFINE MODEL PARAMETERS
% Used MCS
MCS = wifi_params.MCS(i_mcs);

cbpb_index = 2; % normal length, TEST

noiseVar = 10^(-SNR_dB/20); % noise variance, without estimation, perfect kowledge



% Channel coding
% type of PHY layer specific
% auxiliary
switch MCS.phy_type
    case {'Ctrl'}
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb);
        % number of symbol padding bits
        N_blk_pad = (N_blks*MCS.N_cbpb)-(N_cw*L_ldpc);
        
    case {'SC'}
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
        % number of symbol padding bits
        N_blk_pad = (N_blks*MCS.N_cbpb(cbpb_index))-(N_cw*L_ldpc);
        % LDPC decoder definition
        hDec = comm.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied');
        % hDec = comm.gpu.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied');
    case 'OFDM'
        N_sym = wifi_params.coding.N_sym;
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % num of data bits
        num_data_bits = wifi_params.coding.num_data_bits;
        % LDPC decoder definition
        hDec = comm.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied');
        % hDec = comm.gpu.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied');
        
    otherwise
        error('Wrong PHY type.')
end




% Spreading
if strcmp(wifi_params.general.PHYlayer, 'OFDM') == 0
    n_fft = wifi_params.mapping.n_fft; % chips
    n_tot = wifi_params.mapping.n_tot; % chips
    n_gi = n_fft-n_tot; % chips
else
    n_fft = wifi_params.mapping.n_fft; % samples
    n_cp = wifi_params.cprefix.n_cp; % samples
    n_gi = wifi_params.mapping.N_GI; % samples
end

%% RECEIVED VECTOR
% y_rx = ChanMod_output;
y_rx = chanObj.output.x_hn;
%% FRAME PARSER
% extract preambles and other fields from frame
% no synchronization is performed, yet
% load exact number of chips for each frame field (excepting Data and TRN -> here set to -1)
y_rxParsed = parseFrame(y_rx, wifi_params, MCS, chanObj);

%% STF - SYNCHRONIZATION
% using STF
y_rx_stf = y_rxParsed.stf;
% TBD
%% CEF - CHANNEL ESTIMATION
% using CEF
y_rx_cef = y_rxParsed.cef;
% [hEst, HEst] = ad_cir_estimator(y_rx_cef, y_rx_stf, chanObj, wifi_params, MCS);

%% DATA - GUARD REMOVAL,  SYMBOL DE-BLOCKING
% y_rx_data = y_rxParsed.data;
[codedOut, uncodedOut] = ayRxDataField(y_rxParsed, chanObj, wifi_params, MCS, noiseVar);



%% Save results

rxObj.output.PSDU = codedOut.PDSU_rx; % user traffic data without overhead
rxObj.output.PSDU_padded = codedOut.PDSU_rx_padded; % user traffic data without overhead padded with zeros

if strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    rxObj.output.coded_data = codedOut.decoderOut; % decoder output - user traffic data + SERVICE zeros field with overhead
    rxObj.output.uncoded_data = (uncodedOut.dataDemodulatorOutLLR(1:length(uncodedOut.dataDemodulatorOutLLR)-N_blk_pad)<0) + 0; % decoder input (raw data) (demodulated bits)
    rxObj.output.uncoded_dataLLR = uncodedOut.dataDemodulatorOutLLR;
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC ---------------------------------------------------------
    error;
    rxObj.output.coded_data = rs_decoderOut.'; % decoder output - user traffic data + SERVICE zeros field with overhead
    rxObj.output.uncoded_data = deintrlIn.'; % decoder input (raw data) (demodulated bits)
elseif strcmp(MCS.phy_type,'OFDM') == 1
    error;
    rxObj.output.coded_data = decoderOut_padded.'; % decoder output - same as decoderOut in this case
    rxObj.output.uncoded_data = demodulatorOutHD_OFDM + 0; % decoder input (raw data) (demodulated bits)
else
    rxObj.output.coded_data = decoderOut_padded.'; % PSDU_rx_padded; % decoder output
    rxObj.output.uncoded_data = dataDemodulatorOutHD + 0; % decoder input (raw data) (demodulated bits)
end
end