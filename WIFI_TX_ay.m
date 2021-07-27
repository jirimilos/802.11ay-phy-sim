function txObj = WIFI_TX_ay(wifi_params,txObj,mcs_num)
% WIFI transmitter function according to IEEE 802.11ay standard
% Author:   Jiri Milos, DREL FEEC BUT, 2020

%% DEFINE MODEL PARAMETERS
% Used MCS
MCS = wifi_params.MCS(mcs_num);
    

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

%% GENERATE PSDU INFORMATION
PSDU_tx = generatePSDU(wifi_params);

%% PSDU PADDING
[PSDU_tx_padded, PSDU_tx] = paddingPSDU(PSDU_tx, wifi_params, MCS);

%% SCRAMBLING
[scramblerOut, scramblingSeq_data, scramblingSeq_block_pad] = txScrambling(PSDU_tx_padded, wifi_params, MCS);
% save('scramblerOut.mat','scramblerOut')
%% LDPC CHANNEL CODING
[encoderOut, encoderOut_single_row] = ayLDPCEncoding(scramblerOut, wifi_params, MCS, scramblingSeq_block_pad);
% save('encoderOut.mat','encoderOut')
%% MODULATION MAPPING
[modulatedDataSymbol_rotated] = ayModulator(encoderOut, wifi_params, MCS);
% konstelace(modulatedDataSymbol, modulatedDataSymbol_rotated)

%% SYMBOL BLOCKING AND GUARD INSERTION incl. STBC (data part only) + OTHER FIELDS PROCESSING
[x_s, STF, CEF, Header, BTF] = ayBlockingAndGI(modulatedDataSymbol_rotated, wifi_params, MCS);


%% Save outputs
% PSDU
txObj.output.PSDU = PSDU_tx; % user traffic data without overhead
txObj.output.PSDU_padded = PSDU_tx_padded; % user traffic data without overhead padded with zeros
% coded BER
txObj.output.coded_data = scramblerOut; % encoder input = scrambler out
% uncoded (raw BER)
txObj.output.uncoded_data = encoderOut_single_row.'; % encoder output (raw data)

% bit message of other fields
txObj.output.STF = STF;
txObj.output.CEF = CEF;
txObj.output.Header = Header;
txObj.output.BTF = BTF;

% trasmitter output symbols
txObj.output.x_s = x_s;

% necessary genie data
txObj.genie.scramblingSeq_data = scramblingSeq_data;

end