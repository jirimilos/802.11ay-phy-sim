function [PSDU_rx, PSDU_rx_padded, Header_rx] = extractPSDU(descramblerOut, wifi_params, MCS)

% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;

if (strcmp(MCS.phy_type,'SC') == 1) || (strcmp(MCS.phy_type,'OFDM') == 1) % SC or OFDM ---------------------------------------------------------
    PSDU_rx_padded = descramblerOut;
    PSDU_rx = PSDU_rx_padded(1:LENGTH*8);
    Header_rx = [];
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC ---------------------------------------------------------
    PSDU_rx_padded = descramblerOut.';
    PSDU_rx = descramblerOut.';
    Header_rx = [];
else % Ctrl -----------------------------------------------------------
%     PSDU_rx_padded = decoderOut_padded;
    PSDU_rx_padded = descramblerOut;
    PSDU_rx = descramblerOut(1,5*8+1:end);
    Header_rx = descramblerOut(1, 1:5*8);
end

end

