function rx_DataSymbols = ayDataDeblocking(y_rx_data, wifi_params, MCS, chanObj)

n_fft = wifi_params.mapping.n_fft; % samples
n_gi = wifi_params.mapping.N_GI; % samples

nRxAnt = chanObj.params.nRX;

y_rx_data_length = length(y_rx_data(:,:,1));

if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
    % equalizer - test
    Q = 32; % number of chips
    n_fft_spr = 2^(nextpow2(Q)+1);
    p = (n_fft_spr-Q)/2;
    %
    %     cfe_length = wifi_params.framing.ChipLengths(2);
    %     slidingWinInput = [y_rx_data(cfe_length-p+1:end,:); zeros(p,1)].';
    %
    %     n_windows = (length(slidingWinInput)-Q)/Q;
    %
    %     rx_DataSymbols = zeros(n_windows,n_fft_spr);
    %     for iw = 1:n_windows
    %         rx_DataSymbols(iw,:) = slidingWinInput(1,(iw-1)*Q+(1:64));
    %     end
    
    % derotation
    n = (0:y_rx_data_length-1).';
    rx_modulatedDataSymbol_derotated = y_rx_data.*repmat(exp(-1j*pi*n/2),1,1,nRxAnt);
    N_blks = y_rx_data_length/Q;
    % reshaping to blocks with 32 columns (32 chips = 1 modulated symbol)
    rx_DataSymbols = reshape(rx_modulatedDataSymbol_derotated, Q, []).';
    
elseif strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    cbpb_index = 2; % normal length, TEST
    % number of codewords
    N_cw = wifi_params.coding.N_cw;
    % single codeword length
    L_ldpc = wifi_params.coding.L_ldpc;
    % number of symbol blocks
    N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
    N_parallelBlks = N_blks/MCS.N_ss;
    % guard removal and data symbols extract
    % throw away the first guard interval in case of SISO, SIMO, MIMO; processing in the 
    % case of TxD is different
    if strcmpi(wifi_params.antConfig.mode,'TxD') == 1
%         if strcmpi(wifi_params.channel.type,'awgn') == 1
%             y_rx_data_wGI = y_rx_data(n_gi+1:end,:,:);
%         else
            y_rx_data_wGI = y_rx_data(n_gi/2+1:end-(n_gi/2),:,:); % remove the first and the last n_gi/2 symbols
%         end
    else
        y_rx_data_wGI = y_rx_data(n_gi+1:end,:,:);
    end
    
    rx_DataSymbols = zeros(N_parallelBlks, n_fft, nRxAnt); % preallocation
    % reshape per Rx Antenna
    for iRxAnt = 1:nRxAnt
        
        y_rx_data_wGIrows = reshape(y_rx_data_wGI(:,:,iRxAnt),[],N_parallelBlks).'; % streams
        rx_DataSymbols(:,:,iRxAnt) = y_rx_data_wGIrows;
        
    end
    
    
elseif strcmp(MCS.phy_type,'OFDM') == 1 % OFDM ---------------------------------------------------------------
    N_tones = wifi_params.mapping.n_tot.nonzero; % number of active subcarriers
    y_rxOFDM_cp_ifft = reshape(y_rx_data, n_fft+n_cp, []);
    % remove cyclic prefix
    y_rxOFDM_ifft = y_rxOFDM_cp_ifft(n_cp+1:end,:);
    % perform FFT
    y_rx_OFDM_data_pilots = (N_tones/n_fft)*fft(y_rxOFDM_ifft, n_fft); % see Sec. 20.5.3.2.6
    %     y_rx_OFDM_data_pilots = (1/sqrt(N_tones))*fft(y_rxOFDM_ifft, n_fft); % see Sec. 20.5.3.2.6
    % get data symbols
    rx_modulatedDataSymbol_blocks = y_rx_OFDM_data_pilots(wifi_params.mapping.ir_data,:);
    
end

end

