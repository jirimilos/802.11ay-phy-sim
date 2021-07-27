function [scramblerOut, scramblingSeq_data, scramblingSeq_block_pad] = txScrambling(PSDU_tx_padded, wifi_params, MCS)

if (strcmp(MCS.phy_type, 'Ctrl') == 1)
    scramblerOut = PSDU_tx_padded;
    scramblingSeq_data = [];
    scramblingSeq_block_pad = [];
    
elseif (strcmp(MCS.phy_type, 'SC') == 1)
    % number of symbol padding bits
    N_blk_pad = wifi_params.coding.N_blk_pad;
        
    scr_seed = wifi_params.scrambling.scr_seed;
    scramblingSeq = AD_ScramblingSeqGen(length(PSDU_tx_padded)+N_blk_pad, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length(PSDU_tx_padded)); % scrambling seq for input data
    scramblingSeq_block_pad = scramblingSeq(length(PSDU_tx_padded)+1:end); % scrambling sequence for blocks pad bits, see ch. 20.6.3.2.3.3/e, page 2475
    if wifi_params.general.useScrambling == 0
        scramblerOut = xor(PSDU_tx_padded, zeros(size(scramblingSeq_data))); % tx scrambling OFF
    else
        scramblerOut = xor(PSDU_tx_padded, scramblingSeq_data); % tx scrambling ON
    end
    
elseif (strcmp(MCS.phy_type, 'OFDM') == 1)
    scr_seed = wifi_params.scrambling.scr_seed;
    scramblingSeq = AD_ScramblingSeqGen(length(PSDU_tx_padded)+504, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length(PSDU_tx_padded)); % scrambling seq for input data
    if wifi_params.general.useScrambling == 0
        scramblerOut = xor(PSDU_tx_padded, zeros(size(scramblingSeq_data))); % tx scrambling OFF
    else
        scramblerOut = xor(PSDU_tx_padded, scramblingSeq_data); % tx scrambling ON
    end
    
end

end

