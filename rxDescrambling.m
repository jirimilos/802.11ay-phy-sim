function [descramblerOut] = rxDescrambling(descramblerIn, wifi_params, MCS)

cbpb_index = 2; % normal length, TEST

% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;

if strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    % number of codewords
    N_cw = wifi_params.coding.N_cw;
    % single codeword length
    L_ldpc = wifi_params.coding.L_ldpc;
    % number of symbol blocks
    N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
    
    % descrambling sequence defined in line 100
    % generate scrambling sequence of data and block pad length
    scr_seed = wifi_params.scrambling.scr_seed;
    length_PSDU_rx_padded = 8*LENGTH+wifi_params.coding.N_Data_pad;
    % number of symbol padding bits
    N_blk_pad = (N_blks*MCS.N_cbpb(cbpb_index))-(N_cw*L_ldpc);
    
    scramblingSeq = AD_ScramblingSeqGen(length_PSDU_rx_padded+N_blk_pad, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length_PSDU_rx_padded); % scrambling seq for input data
    
    if wifi_params.general.useScrambling == 0
        descramblerOut = xor(descramblerIn, zeros(size(scramblingSeq_data.')));
    else
        descramblerOut = xor(descramblerIn, scramblingSeq_data);
    end
    
elseif strcmp(MCS.phy_type, 'OFDM') == 1 % OFDM --------------------------------------------------------
    
    scr_seed = wifi_params.scrambling.scr_seed;
    scramblingSeq = AD_ScramblingSeqGen(length(descramblerIn)+504, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length(descramblerIn)).'; % scrambling seq for input data
    if wifi_params.general.useScrambling == 0
        descramblerOut = xor(descramblerIn, zeros(size(scramblingSeq_data))); % tx scrambling OFF
    else
        descramblerOut = xor(descramblerIn, scramblingSeq_data); % tx scrambling ON
    end
    
    % elseif strcmp(MCS.phy_type, 'LPSC') == 1 % LPSC --------------------------------------------------------
    %     scr_seed = wifi_params.scrambling.scr_seed;
    %     descramblingSeq = AD_ScramblingSeqGen(N_eb+N_blk_pad, scr_seed);
    %     descramblingSeq_data = descramblingSeq(1:length(rs_decoderOut)); % descrambling seq for input data
    %
    %     if wifi_params.general.useScrambling == 0
    %         descramblerOut = rs_decoderOut; % tx scrambling OFF
    %     else
    %         descramblerOut = xor(rs_decoderOut, descramblingSeq_data); % tx scrambling ON
    %     end
    
else % Ctrl -----------------------------------------------------------
    descramblerOut = descramblerIn;
end


end

