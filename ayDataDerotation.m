function rx_DataSymbols = ayDataDerotation(y_rx_data, wifi_params, MCS, chanObj)

n_fft = wifi_params.mapping.n_fft; % samples
n_gi = wifi_params.mapping.N_GI; % samples

nRxAnt = chanObj.params.nRX;

y_rx_data_length = length(y_rx_data(:,:,1));
if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
    % derotation
    n = (0:y_rx_data_length-1).';
    rx_modulatedDataSymbol_derotated = y_rx_data.*repmat(exp(-1j*pi*n/2),1,1,nRxAnt);
    N_blks = y_rx_data_length/32;
    % reshaping to blocks with 32 columns (32 chips = 1 modulated symbol)
    rx_DataSymbols = reshape(rx_modulatedDataSymbol_derotated, 32, []).';
    %
elseif strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    cbpb_index = 2; % normal length, TEST
    % number of codewords
    N_cw = wifi_params.coding.N_cw;
    % single codeword length
    L_ldpc = wifi_params.coding.L_ldpc;
    % number of symbol blocks
    N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
        
    % guard removal and data symbols extract
    % throw away the first guard interval
    y_rx_data_wGI = y_rx_data(n_gi+1:end,:);
    
    
    %%%% UPRAVA 31/03/2020:
    % odstraneni derotace -> bude provedena az po ekvalizaci!!!
    
    y_rx_data_wGIrows = reshape(y_rx_data_wGI,[],N_blks).';
    rx_DataSymbols = y_rx_data_wGIrows;
    
%     
%     % data part
%     modulated_symbol_index_k = (0:n_fft-n_gi-1);
% %     modulated_symbol_derotating_factor = repmat(exp(-1j*pi*modulated_symbol_index_k/2),1,1,nRxAnt);
%     dataPartDerotMatrix = repmat(exp(-1j*pi*modulated_symbol_index_k/2),N_blks,1);
%     
%     % GI part
%     k_Ga64 = 0:n_gi-1;
% %     repmat((repmat(exp(-1j*pi*k_Ga64/2), 1, N_blks+1)).',1,1,nRxAnt);
%     giPartDerotMatrix = repmat(exp(-1j*pi*k_Ga64/2),N_blks,1);
%     
%     % cancatenate derotation matrix
%     derotMatrix = [dataPartDerotMatrix, giPartDerotMatrix];
%     if wifi_params.general.dataSymbolRotation == 0
%         rx_DataSymbols = y_rx_data_wGIrows;
%     else
%         rx_DataSymbols = y_rx_data_wGIrows.*derotMatrix; % deratoted data+gi symbols [N_blks x 512]
%     end
    
    
% %     rx_DataSymbols_rotated = zeros((n_fft-n_gi)*N_blks,1,nRxAnt); % preallocation (N_blks of Data symbols within data frame) 
% %     rx_Ga64_rotated = zeros(n_gi*(N_blks+1),1,nRxAnt); % preallocation (N_blks +1 of Ga64 sequences within data frame)
% %     
% %     idxData = kron((64:n_fft:length(y_rx_data(:,:,1))-1), ones(1,448))+repmat(1:448,1,N_blks); % GI indices per Rx antenna
% %     idxGI = kron((1:n_fft:length(y_rx_data(:,:,1)))-1, ones(1,64))+repmat(1:64,1,N_blks+1); % GI indices per Rx antenna
% %     
% %     rx_DataSymbols_rotated(:,1,:) = y_rx_data(idxData,1,:);
% %     rx_Ga64_rotated(:,1,:) = y_rx_data(idxGI,1,:);
% %     
% %     if mod(length(y_rx_data(:,:,1))-64,n_fft) ~= 0
% %         error('Wrong length of block of received symbols');
% %     end
%     
% %     rx_modulatedDataSymbolsAndGuard = reshape(y_rx_data, n_fft, []).';
% %     % rest of guard removal
% %     rx_guard = rx_modulatedDataSymbolsAndGuard(:,1:n_gi);
% %     rx_Ga64_rotated(1:N_blks, :) = rx_guard;
% %     k_Ga64 = 0:n_gi-1;
% %     rx_Ga64 = rx_Ga64_rotated.*repmat((repmat(exp(-1j*pi*k_Ga64/2), 1, N_blks+1)).',1,1,nRxAnt);
%     
% %     rx_modulatedDataSymbol_blocks = rx_modulatedDataSymbolsAndGuard(:,n_gi+1:n_fft);
%     
%     % derotate data symbols first
% %     modulated_symbol_index_k = (0:numel(rx_DataSymbols_rotated(:,1,1))-1).';
% %     modulated_symbol_derotating_factor = repmat(exp(-1j*pi*modulated_symbol_index_k/2),1,1,nRxAnt);
% %     if wifi_params.general.dataSymbolRotation == 0
% %         rx_DataSymbols = rx_DataSymbols_rotated;
% %     else
% %         rx_DataSymbols = rx_DataSymbols_rotated.*modulated_symbol_derotating_factor;
% %     end
    
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

