function y_rxParsed = parseFrame(y_rx, wifi_params, MCS, chanObj)

nRxAnt = chanObj.params.nRX;
% parsing
switch MCS.phy_type
    case 'Ctrl' % Control -----------------------------------------------------------------
        if nRxAnt > 1
            error('Only single receiving antenna for Ctrl physical layer');
        end
        y_rx_stf    = y_rx(wifi_params.framing.ChipMap{1});
        y_rx_cef    = y_rx(wifi_params.framing.ChipMap{2});
        y_rx_header = [];
        y_rx_data   = y_rx(wifi_params.framing.ChipMap{3}); % header + data
        y_rx_btf    = y_rx(wifi_params.framing.ChipMap{4});
    case {'SC', 'LPSC', 'OFDM'} % SC, LPSC or OFDM ------------------------------------------------------
        if wifi_params.general.sendDataOnly == 0
            y_rx_stf    = y_rx(find(wifi_params.framing.ChipMap{1} == 1),:,:);
            y_rx_cef    = y_rx(find(wifi_params.framing.ChipMap{2} == 1),:,:);
            y_rx_header = y_rx(find(wifi_params.framing.ChipMap{3} == 1),:,:);
            y_rx_data   = y_rx(find(wifi_params.framing.ChipMap{4} == 1),:,:);
            y_rx_btf    = y_rx(find(wifi_params.framing.ChipMap{5} == 1),:,:);
        else
            y_rx_stf    = [];
            y_rx_cef    = [];
            y_rx_header = [];
            y_rx_data   = y_rx;
            y_rx_btf    = [];
        end
    otherwise
        error('Not defined yet');
end

y_rxParsed.stf = y_rx_stf;
y_rxParsed.cef = y_rx_cef;
y_rxParsed.header = y_rx_header;
y_rxParsed.data = y_rx_data;
y_rxParsed.btf = y_rx_btf;

end

