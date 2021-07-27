function [PSDU_tx_padded, PSDU_tx] = paddingPSDU(PSDU_tx, wifi_params, MCS)

% number of codewords
N_cw = wifi_params.coding.N_cw;

if (strcmp(MCS.phy_type,'SC') == 1) || (strcmp(MCS.phy_type,'OFDM') == 1)
    padding_zeros = zeros(wifi_params.coding.N_Data_pad,1);
    PSDU_tx_padded = [PSDU_tx; padding_zeros];

elseif strcmp(MCS.phy_type,'LPSC') == 1
    PSDU_tx_padded = PSDU_tx;

elseif strcmp(MCS.phy_type,'Ctrl') == 1
    % Header is encoded together with the user data, see 20.4.3.2.3
    % Generate Header bits
    % (provisionally) generated as a true random data
    Header = randi([0 1], 5*8, 1); % 5 octets as in Table 20-11, page 2454
    HeaderCRC = Header(25:end);
    
    PSDU_tx = cell(1, N_cw);
    PSDU_tx_padded = cell(1, N_cw);
    for i_cw = 1:N_cw
        if i_cw == 1
            PSDU_tx_tmp = randi([0 1], 1, wifi_params.coding.num_data_bits(i_cw)-length(Header));
            PSDU_tx_tmp_padded = [Header.', PSDU_tx_tmp, zeros(1, wifi_params.coding.N_Data_pad(i_cw))];
        else
            PSDU_tx_tmp = randi([0 1], 1, wifi_params.coding.num_data_bits(i_cw));
            PSDU_tx_tmp_padded = [PSDU_tx_tmp, zeros(1, wifi_params.coding.N_Data_pad(i_cw))];
        end
        PSDU_tx{i_cw} = PSDU_tx_tmp;
        PSDU_tx_padded{i_cw} = PSDU_tx_tmp_padded;
    end
    PSDU_tx = [PSDU_tx{:}].';
    PSDU_tx_padded = [PSDU_tx_padded{:}].';
    
    if numel(PSDU_tx_padded) ~= 168*N_cw
        error('Error in Header + Data length definition!');
    end
    
end
% Ctrl: Scrambling defined in Sec. 20.4.3.3.2, page 2454
% SC: both should be scrambled, padding zeros with continuous scrambling word, see page 2474

end

