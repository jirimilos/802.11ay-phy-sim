function [codedWordStruct, uncodedWordStruct] = ayRxDataField(y_rxParsed, chanObj, wifi_params, MCS, noiseVar)

switch MCS.phy_type
    case {'Ctrl'} % Control -----------------------------------------------------------------
        y_rx_data = [y_rxParsed.cef; y_rxParsed.data]; % in this case we need last p symbols of CEF part (overlap-cut method)
    case {'SC'}
        y_rx_data = y_rxParsed.data;
    otherwise
        
end


% data symbols derotation
rx_modulatedDataAndGISymbol_blocks = ayDataDeblocking(y_rx_data, wifi_params, MCS, chanObj);

%% CHANNEL ESTIMATION
% perfect channel knowledge or estimated channel
if strcmp(wifi_params.general.channelKnowledge,'perfect') == 1
    % get perfect CTF
    HEst = chanObj.genie.H;
    hEst = chanObj.genie.h;
end

%% DATA BLOCKS DEMODULATION
dataDemodulatorOutLLR = ayDemodulator(rx_modulatedDataAndGISymbol_blocks, chanObj, HEst, wifi_params, MCS, noiseVar);
dataDemodulatorOutHD = (dataDemodulatorOutLLR < 0);


%% CHANNEL DECODING
[decoderOut] = ayLDPCDecoding(dataDemodulatorOutLLR, wifi_params, MCS);

%% DESCRAMBLING
descramblerOut = rxDescrambling(decoderOut, wifi_params, MCS);

%% EXTRACT PSDU INFORMATION
% PSDU, PSDU padded and Header_rx (for CTRL only)
[PSDU_rx, PSDU_rx_padded, Header_rx] = extractPSDU(descramblerOut, wifi_params, MCS);

codedWordStruct.PDSU_rx = PSDU_rx;
codedWordStruct.PDSU_rx_padded = PSDU_rx_padded;
codedWordStruct.Header_rx = Header_rx;

codedWordStruct.decoderOut = decoderOut;
codedWordStruct.descramblerOut = descramblerOut;

uncodedWordStruct.dataDemodulatorOutHD = dataDemodulatorOutHD;
uncodedWordStruct.dataDemodulatorOutLLR = dataDemodulatorOutLLR;
end

