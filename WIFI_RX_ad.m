function rxObj = WIFI_RX_ad(wifi_params,rxObj,channelObj, ChanMod_output, mcs_num, SNR_dB, txObj, var_n, var_s)
% WIFI transmitter function according to IEEE 802.11ad standard
% Author:   Jiri Milos, DREL FEEC BUT, 2019

%% DEFINE MODEL PARAMETERS
% Used MCS
MCS = wifi_params.MCS(mcs_num);


% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;
% Channel coding
% type of PHY layer specific
switch MCS.phy_type
    case {'Ctrl', 'SC'}
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb);
        % number of symbol padding bits
        N_blk_pad = (N_blks*MCS.N_cbpb)-(N_cw*L_ldpc);
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
    case 'LPSC'
        % inner Reed Solomon
        % number of codewords
        N_cw = wifi_params.coding.RS.N_cw;
        N_encsym_tot = wifi_params.coding.RS.N_encsym_tot;
        n_rs_encoders = wifi_params.coding.RS.numEncoders;
        hEnc_RS = cell(1, n_rs_encoders);
        % RS encoder(s) definition
        for ienc = 1:n_rs_encoders
            n = wifi_params.coding.RS.n_k(ienc, 1);
            k = wifi_params.coding.RS.n_k(ienc, 2);
            hDec_RS{ienc} = comm.RSDecoder(...
                'BitInput', false,...
                'CodewordLength', n,...
                'MessageLength', k,...
                'PrimitivePolynomialSource', 'Property',...
                'PrimitivePolynomial', wifi_params.coding.RS.P_bin,...
                'GeneratorPolynomialSource','Property',...
                'GeneratorPolynomial', wifi_params.coding.RS.G{ienc});
        end
        
        % outer Block
        % number of symbol blocks
        N_blks = wifi_params.coding.Blck.N_blks;
        % number of symbol padding bits
        N_blk_pad = wifi_params.coding.Blck.N_blk_pad;
        blck_genmat = wifi_params.coding.Blck.G;
        blck_n = wifi_params.coding.Blck.n_k(1);
        blck_k = wifi_params.coding.Blck.n_k(2);
        N_eb = wifi_params.coding.Blck.N_eb;
        genmat = wifi_params.coding.Blck.G;
    otherwise
        error('Wrong PHY type.')
end

% Modulation mapping
modNorm = wifi_params.modulation(log2(MCS.M)).hMod_factor;
hDemod = wifi_params.modulation(log2(MCS.M)).hDemod;
hDemod_LLR = wifi_params.modulation(log2(MCS.M)).hDemod_LLR;

SymbolAlphabet = (constellation(hDemod).')*wifi_params.modulation(log2(MCS.M)).hMod_factor;
hScalQuant = wifi_params.modulation.hScalQuant;

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


Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
%% RECEIVED VECTOR
y_rx = ChanMod_output;

%% FRAME PARSER
% extract preambles and other fields from frame
% no synchronization is performed, yet
% load exact number of chips for each frame field (excepting Data and TRN -> here set to -1)

% Control

% SC

% LPSC

% OFDM


% parsing
switch MCS.phy_type
    case 'Ctrl' % Control -----------------------------------------------------------------
        y_rx_stf = y_rx(wifi_params.framing.ChipMap{1});
        y_rx_cef = y_rx(wifi_params.framing.ChipMap{2});
        y_rx_data = y_rx(wifi_params.framing.ChipMap{3}); % header + data
        y_rx_btf = y_rx(wifi_params.framing.ChipMap{4});
    case {'SC', 'LPSC', 'OFDM'} % SC, LPSC or OFDM ------------------------------------------------------
        y_rx_stf = y_rx(wifi_params.framing.ChipMap{1});
        y_rx_cef = y_rx(wifi_params.framing.ChipMap{2});
        y_rx_header = y_rx(wifi_params.framing.ChipMap{3});
        y_rx_data = y_rx(wifi_params.framing.ChipMap{4});
        y_rx_btf = y_rx(wifi_params.framing.ChipMap{5});
    otherwise
        error('Not defined yet');
end

%% SYNCHRONIZATION
% using STF
% TBD

%% CHANNEL ESTIMATION
% using CEF
est_h = ad_cir_estimator(y_rx_cef, y_rx_stf, wifi_params);

%% DATA - GUARD REMOVAL AND SYMBOL DE-BLOCKING
y_rx_data_length = length(y_rx_data);
if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
    % derotation
    n = 0:y_rx_data_length-1;
    rx_modulatedDataSymbol_derotated = y_rx_data.*exp(-1j*pi*n/2);
    N_blks = y_rx_data_length/32;
    % reshaping to blocks with 32 columns (32 chips = 1 modulated symbol)
    rx_modulatedDataSymbol_blocks_derotated = reshape(rx_modulatedDataSymbol_derotated, 32, []).';
    %
elseif strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    % guard removal
    rx_Ga64_rotated = zeros(N_blks+1, n_gi); % preallocation (N_blks +1 of Ga64 sequences within data frame)
    % remove Ga64 at the end of rx data blocks stream
    
    rx_Ga64_rotated(N_blks+1, :) = y_rx_data(1, y_rx_data_length-n_gi+1:y_rx_data_length);
    y_rx_data(y_rx_data_length-n_gi+1:y_rx_data_length) = [];
    
    if mod(length(y_rx_data),n_fft) ~= 0
        error('Wrong length of block of received symbols');
    end
    
    rx_modulatedDataSymbolsAndGuard = reshape(y_rx_data, n_fft, []).';
    % rest of guard removal
    rx_guard = rx_modulatedDataSymbolsAndGuard(:,1:n_gi);
    rx_Ga64_rotated(1:N_blks, :) = rx_guard;
    k_Ga64 = 0:n_gi-1;
    rx_Ga64 = rx_Ga64_rotated.*repmat(exp(-1j*pi*k_Ga64/2), N_blks+1, 1);
    
    rx_modulatedDataSymbol_blocks = rx_modulatedDataSymbolsAndGuard(:,n_gi+1:n_fft);
    % derotate data symbols first
    modulated_symbol_index_k = (0:numel(rx_modulatedDataSymbol_blocks(1,:))-1);
    modulated_symbol_derotating_factor = exp(-1j*pi*modulated_symbol_index_k/2);
    if wifi_params.general.dataSymbolRotation == 0
        rx_modulatedDataSymbol_blocks_derotated = rx_modulatedDataSymbol_blocks;
    else
        rx_modulatedDataSymbol_blocks_derotated = rx_modulatedDataSymbol_blocks.*repmat(modulated_symbol_derotating_factor, N_blks, 1);
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
    
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC ---------------------------------------------------------------
    % guard removal (the last block )a64
    y_rx_data(end-64+1:end) = []; % erase the last Ga64 block (get the field lenght as a multiple of 512)
    % create subblocks of 64 symbols
    y_rx_data_subblocks = reshape(y_rx_data,64,[]).';
    % get d56blocks and G8
    y_rx_data_G64 = y_rx_data_subblocks((((1:wifi_params.coding.Blck.N_blks)-1)*8)+1,:);
    y_rx_data_subblocks((((1:wifi_params.coding.Blck.N_blks)-1)*8)+1,:) = []; % erase G64 blocks
    y_rx_data_d56andG8 = y_rx_data_subblocks;
    % pick d56 blocks (G8 unused)
    d56_blocks = y_rx_data_d56andG8(:,1:56);
    d56_blocks_tmp = d56_blocks.';
    rx_modulatedDataSymbol_rotated = d56_blocks_tmp(:).';
    modulated_symbol_index_k = (0:numel(rx_modulatedDataSymbol_rotated(1,:))-1);
    modulated_symbol_derotating_factor = exp(-1j*pi*modulated_symbol_index_k/2);
    % data symbols derotation
    rx_modulatedDataSymbol_blocks = rx_modulatedDataSymbol_rotated.*modulated_symbol_derotating_factor;
end
%% CHANNEL EQUALIZATION
% with LS ZF equalizer
% get channel CIR and CFT
switch MCS.phy_type
    case 'Ctrl' % Control -----------------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            est_h_padded = [channelObj.h, zeros(1, 32-length(channelObj.h))];
            H_est = channelObj.H(1:32);
        else
            est_h_padded = est_h(1:32);
            H_est = fft(est_h_padded);
        end
        
    case 'SC' % SC ------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            est_h_padded = [channelObj.h, zeros(1, 448-length(channelObj.h))];
            H_est = channelObj.H(1:end-n_gi);
        else
            est_h_padded = [[est_h(1:120),zeros(1,length(est_h)-120)], zeros(1,n_fft-n_gi-length(est_h))];
            H_est = fft(est_h_padded);
        end
        
    case 'OFDM' % OFDM ------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            est_h_padded = [channelObj.h, zeros(1, n_fft-length(channelObj.h))];
            H_est = channelObj.H(wifi_params.mapping.ir_data,:); % 336 symbol positions
        else
            est_h_padded = [[est_h(1:120),zeros(1,length(est_h)-120)], zeros(1,n_fft-length(est_h))];
%             H_est = fft(est_h_padded);
%             H_est = H_est(wifi_params.mapping.ir_data);
%             H_est = H_est(1:336);
            H_est = channelObj.H(wifi_params.mapping.ir_data);
        end
        
    case 'LPSC'
        %         H_est % ma byt bez 8*7 symbolu (mozna bude potreba i jina uprava)
        if strcmp(wifi_params.general.channel,'awgn') == 1
            est_h_padded = [channelObj.h, zeros(1, 56*7-length(channelObj.h))];
            H_est = channelObj.H(1:56*7);
        else
            est_h_padded = [[est_h(1:120),zeros(1,length(est_h)-120)], zeros(1,56*7-128)];
            H_est = fft(est_h_padded);
        end
    otherwise % OFDM -----------------------------------------------------------------
        
end

% get perfect CTF
H_est_perf = channelObj.H;

% other preprocessing
if strcmp(MCS.phy_type,'Ctrl') == 1 % Control ------------------------------------------------------------
    FDE_data_output = zeros(size(rx_modulatedDataSymbol_blocks_derotated));
    rx_layer_x = zeros(N_blks, 1);
    rx_user_symbols = zeros(N_blks, 1);
    demodulatorOutSSD = zeros(N_blks, log2(MCS.M)*length(FDE_data_output(1,:)));
    
elseif strcmp(MCS.phy_type,'SC') == 1 % SC ---------------------------------------------------------------
    FDE_data_output = zeros(size(rx_modulatedDataSymbol_blocks));
    rx_layer_x = zeros(N_blks, length(FDE_data_output(1,:)));
    demodulatorOutSSD = zeros(N_blks, log2(MCS.M)*length(FDE_data_output(1,:)));
    
elseif strcmp(MCS.phy_type,'OFDM') == 1 % OFDM ---------------------------------------------------------------
    FDE_data_output = zeros(size(rx_modulatedDataSymbol_blocks));
    rx_layer_x = zeros(length(FDE_data_output(:,1)),N_sym);
    demodulatorOutSSD = zeros(N_sym, log2(MCS.M)*length(FDE_data_output(1,:)));
    
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC -----------------------------------------------------------
    FDE_data_output = zeros(size(rx_modulatedDataSymbol_blocks));
    rx_layer_x = zeros(1, length(FDE_data_output(1,:)));
    demodulatorOutSSD = zeros(1, log2(MCS.M)*length(FDE_data_output(1,:)));
    
else
    % error();
end


% processing per blocks in SC, OFDM and Ctrl, and in signle row for LPSC
switch MCS.phy_type
    case {'Ctrl', 'SC'}
        for i_eq = 1:N_blks
            
            if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
                % prealloc
                FDE_data_input = (rx_modulatedDataSymbol_blocks_derotated(i_eq,:));
                rx_layer_x_tmp = zeros(1, 32);
                rx_user_symbols_tmp = zeros(1, 32);
                for jj = 1:length(FDE_data_input)
                    rx_user_sym_tmp = fft(FDE_data_input(jj));
                    H_est_tmp = H_est(jj);
                    
                    rx_layer_x_tmp(1, jj) = (ifft((H_est_tmp.^(-1)).*rx_user_sym_tmp))*Ga32(jj); % ZF + despreading using Ga32
                    rx_user_symbols_tmp(1, jj) = FDE_data_input(jj)*Ga32(jj); % for SSD
                end
                rx_layer_x(i_eq, 1) = sum(rx_layer_x_tmp)/32;
                rx_user_symbols(i_eq, 1) = sum(rx_user_symbols_tmp)/32;
                H_est_tmp2 = sum(H_est)/32;
            elseif strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
                % prealloc
                demodulatorOutSSD_blk = zeros(1, log2(MCS.M)*length(FDE_data_output(1,:)));
                %
                CSM = hDemod.CustomSymbolMapping;
                if MCS.M == 2
                    bittable = (logical(de2bi(0:length(SymbolAlphabet)-1).'));
                else
                    bittable = flipud(logical(de2bi(CSM).'));
                end
                
                FDE_data_input = (rx_modulatedDataSymbol_blocks_derotated(i_eq,:));
                
                for jj = 1:length(FDE_data_input)
                    
                    rx_user_sym_tmp = fft(FDE_data_input(jj));
                    H_est_tmp = H_est(jj);
                    
                    rx_layer_x(i_eq, jj) = ifft((H_est_tmp.^(-1)).*rx_user_sym_tmp); % ZF
                    %
                    [Q, R] = qr(H_est_tmp);
                    demodulatorOutSSD_blk(1, ((jj-1)*log2(MCS.M))+1:jj*log2(MCS.M)) = (LTE_softsphere(rx_layer_x(i_eq, jj), (rx_user_sym_tmp), Q, R, SymbolAlphabet, bittable, 1, log2(MCS.M))).';
                    
                end
                demodulatorOutSSD(i_eq,:) = demodulatorOutSSD_blk;
                %%% see: http://www.ursi.org/proceedings/procGA11/ursi/C02-3.pdf % jednoducha implementace MMSE
            else % LPSC or OFDM -----------------------------------------------------------------
                
            end
            
        end
        
        if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
            % Differential BPSK
            % prealloc
            demodulatorOutSSD_blk = zeros(1, log2(MCS.M)*length(FDE_data_output(:,1)));
            bittable = (logical(de2bi(0:length(SymbolAlphabet)-1).'));
            for i_k = 1:length(rx_layer_x)
                
                [Q, R] = qr(H_est_tmp2);
                demodulatorOutSSD_blk(1, ((i_k-1)*log2(MCS.M))+1:i_k*log2(MCS.M)) = (LTE_softsphere(rx_layer_x(i_k, 1), rx_user_symbols(i_k, 1), Q, R, SymbolAlphabet, bittable, 1, log2(MCS.M))).';
                
            end
            
            % differential demapping
            
            s_k = zeros(size(demodulatorOutSSD_blk));
            
            for i_k = 1:length(s_k)
                if i_k == 1
                    d_km1 = 1;
                else
                    d_km1 = sign(demodulatorOutSSD_blk(i_k-1));
                end
                s_k( i_k) = demodulatorOutSSD_blk(i_k)*d_km1;
            end
            
            demodulatorOutSSD = s_k;
            %     s_k = 2*encoderOut-1;
            demodulatorOutHD_Ctrl = (demodulatorOutSSD > 0).';
        end
    case 'OFDM'
        demodulatorOutSSD = [];
        for i_eq = 1:N_sym
            
            % bittable and symbol alphabet
            SymbolAlphabet = [1+1j, -1+1j, -1-1j, 1-1j]/sqrt(2);
            bittable = flipud(logical(de2bi([3 1 0 2]).'));
            % samples of whole OFDM  symbol
            FDE_data_input = rx_modulatedDataSymbol_blocks(:,i_eq);
            FDE_data_input_tmp = zeros(size(FDE_data_input));
            % processing within single OFDM symbol
            if (strcmp(MCS.MCS, 'MCS13') == 1) || (strcmp(MCS.MCS, 'MCS14') == 1) % SQPSK is used
                FDE_inp_length = length(FDE_data_input)/2;
                % prealloc
                demodulatorOutSSD_blk = zeros(1, log2(MCS.M)*FDE_inp_length);
                for jj = 1:FDE_inp_length
                    % first and second SQPSK received symbol
                    rx_user_sym_tmp1 = FDE_data_input(jj);
                    rx_user_sym_tmp2 = conj(FDE_data_input(jj+FDE_inp_length));
                    rx_user_sym_tmp = [rx_user_sym_tmp1, rx_user_sym_tmp2];
                    % corresponding H_est
                    H_est_tmp1 = H_est(jj);
                    H_est_tmp2 = H_est(jj+FDE_inp_length);
                    H_est_tmp = [H_est_tmp1, H_est_tmp2]; % complete matrix
                    % first and second SQPSK received symbol - equalized
                    rx_layer_x1 = (H_est_tmp1.^(-1)).*rx_user_sym_tmp1; %
                    rx_layer_x2 = (H_est_tmp2.^(-1)).*rx_user_sym_tmp2; %
                    rx_layer_x = [rx_layer_x1, rx_layer_x2]; % % complete matrix
                    % QR decomposition
                    [Q, R] = qr(H_est_tmp);
                    demodulatorOutSSD_blk(1,((jj-1)*2)+([1, 2])) = (LTE_softsphere(rx_layer_x, rx_user_sym_tmp, Q, R, SymbolAlphabet, bittable, 1, log2(MCS.M)));
                end
                demodulatorOutSSD = [demodulatorOutSSD, demodulatorOutSSD_blk];
            else % QPSK, 16QAM and 64QAM
                FDE_inp_length = length(FDE_data_input);
                kR = FDE_inp_length/2;
                Nsd = 336;
                % prealloc
                demodulatorOutSSD_blk = zeros(1, 2*MCS.N_bpsc*kR);
                %
                %                 for i_c = 1:kR
%                 H_est_tmp = ones(2,1);
                if MCS.M == 4
                    for i_c = 1:kR                        
                        % Static tone pairing demapping
                        d0R = FDE_data_input(i_c);
                        d1R = FDE_data_input(i_c+kR);
                        dR = [d0R; d1R];
                        H_est_tmp = [H_est(i_c); H_est(i_c+kR)];
                        
                        [c_outSSD, ~] = DCM_demod(dR, H_est_tmp, 'QPSK');
                        % save to block
                        demodulatorOutSSD_blk(1, (i_c-1)*4+(1:4)) = c_outSSD;
                    end
                elseif MCS.M == 16
                    for i_c = 1:kR                        
                        % Static tone pairing demapping
                        d0R = FDE_data_input(i_c);
                        d1R = FDE_data_input(i_c+kR);
                        dR = [d0R; d1R];
                        
                        [c_outSSD, ~] = DCM_demod(dR, H_est_tmp(1), '16QAM');
                        % output bits demapping, see 20.5.3.2.4.4
                        % save to block
                        demodulatorOutSSD_blk(1, (i_c-1)*4+(1:4)) = c_outSSD(1:4);
                        demodulatorOutSSD_blk(1, (i_c-1)*4+(1:4)+(2*Nsd)) = c_outSSD(5:8);
                    end
                elseif MCS.M == 64
                    k64_vec = 0:Nsd-1;
                    m64_vec = (112*mod(k64_vec,3))+floor(k64_vec/3); % k mapped to m, as described in Sec. 20.5.3.2.4.5
                    % due to non DCM mapping in 64QAM modulation, the for loop
                    % has doubled number of iterations (for each subcarrier)
                    % interleaving between 3 LDPC codewords on a subcarrier basis
                    % reshape to groups of 6 bits
                    % prealloc
                    c_outSSD = zeros(log2(MCS.M), Nsd);
                    for i_c = 1:Nsd
                        dR = FDE_data_input(i_c);
                        [c_outSSD_tmp, ~] = DCM_demod(dR, H_est_tmp(1), '64QAM');
                        % demapping/deinterleaving bits to correct position
                        c_outSSD(:,m64_vec(i_c)+1) = c_outSSD_tmp;
                    end
                    demodulatorOutSSD_blk = reshape(c_outSSD,1,[]);
                else
                    error('Wrong DCM modulation type');
                end
                
                demodulatorOutSSD = [demodulatorOutSSD, demodulatorOutSSD_blk];
            end
            
        end
        demodulatorOutHD_OFDM = (demodulatorOutSSD > 0);
        
    case 'LPSC'
        FDE_data_input = rx_modulatedDataSymbol_blocks;
        % get modulation bittable
        CSM = hDemod.CustomSymbolMapping;
        if MCS.M == 2
            bittable = (logical(de2bi(0:length(SymbolAlphabet)-1).'));
        else
            bittable = flipud(logical(de2bi(CSM).'));
        end
        
        
        for jj = 1:length(FDE_data_input)
            
            rx_user_sym_tmp = fft(FDE_data_input(jj));
            H_est_tmp = H_est(mod(jj-1,56*7)+1);
            
            rx_layer_x(1, jj) = ifft((H_est_tmp.^(-1)).*rx_user_sym_tmp); % ZF
            %
            [Q, R] = qr(H_est_tmp);
            demodulatorOutSSD(1, ((jj-1)*log2(MCS.M))+1:jj*log2(MCS.M)) = (LTE_softsphere(rx_layer_x(1, jj), (rx_user_sym_tmp), Q, R, SymbolAlphabet, bittable, 1, log2(MCS.M))).';
            
        end
        
    otherwise
        % error()
        
end
% rx_modulatedDataSymbol_blocks_derotated_equalized = FDE_data_output;
rx_modulatedDataSymbol_blocks_derotated_equalized = rx_layer_x;

%% DATA BLOCKS DEMODULATION
demodulatorOut = reshape(-demodulatorOutSSD.',1,[]).';


%% CHANNEL DECODING
decoderIn = 10*demodulatorOut*(1/wifi_params.modulation(log2(MCS.M)).hMod_factor);

if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
    L_cwd = wifi_params.coding.L_cwd;
    N_Data_pad = wifi_params.coding.N_Data_pad;
    
    decoderOut = cell(1, N_cw);
    decoderOut_padded = cell(1, N_cw);
    decoderIn_cw_pattern = zeros(N_cw, L_cwd);
    num_data_bits = wifi_params.coding.num_data_bits;
    for i_cw = 1:N_cw
        
        decoderIn_cw_tmp = decoderIn(1:num_data_bits(i_cw));
        decoderIn(1:num_data_bits(i_cw)) = [];
        decoderIn_cw_pattern_tmp = decoderIn(1:L_cwd);
        decoderIn(1:L_cwd) = [];
        % add zeros to each codeword parity (logical 0 = +1)) for using with
        % LDPC CR = 3/4 (the "strong" logical zero is used; multiplied by 10)
        decoderIn_cw = [decoderIn_cw_tmp; 10*ones(N_Data_pad(i_cw), 1); 10*ones(2*L_cwd, 1); decoderIn_cw_pattern_tmp];
        decoderOut_cw = hDec.step(decoderIn_cw);
        decoderOut{1, i_cw} = decoderOut_cw(1:num_data_bits(i_cw)).';
        decoderOut_padded{1, i_cw} = [decoderOut_cw(1:num_data_bits(i_cw)).', zeros(1, N_Data_pad(i_cw))];
    end
    decoderOut = [decoderOut{:}];
    decoderOut_padded = [decoderOut_padded{:}];
    
elseif strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    % generate scrambling sequence of data and block pad length
    scr_seed = wifi_params.scrambling.scr_seed;
    length_PSDU_rx_padded = 8*LENGTH+wifi_params.coding.N_Data_pad;
    
    scramblingSeq = AD_ScramblingSeqGen(length_PSDU_rx_padded+N_blk_pad, scr_seed);
    
    scramblingSeq_data = scramblingSeq(1:length_PSDU_rx_padded); % scrambling seq for input data
    scramblingSeq_block_pad = scramblingSeq(length_PSDU_rx_padded+1:end); % scrambling sequence for blocks pad bits, see ch. 20.6.3.2.3.3/e, page 2475
    
    
    
    % ee] reverse concatenation of coded bit stream and N_blk_pad zeros
    % demux decoderIn and scrambled block pad zeros (reverse concatenation)
    % extract decoder input bits
    decoderIn_extracted = decoderIn(1:length(decoderIn)-N_blk_pad,:);
    % extract block pad zeros
    extracted_blk_pads = decoderIn(end-N_blk_pad+1:end,:);
    descrambled_blk_pads = extracted_blk_pads.*(-2*scramblingSeq_block_pad+1); % should be all logical zeros (i.e. > 0)
    
    % cc] reshape to matrix with N_blks rows
    % dec_block_matrix = reshape(decoderIn_extracted,[],N_blks).'; % not necessary for decoding
    
    % LDPC decoding (encoding process is fully described in 20.6.3.2.3.3)
    % bb] procedure for converting the scrambled PSDU data to LDPC codeword
    switch MCS.Repetition
        case 1
            if MCS.CR ~= 7/8
                % reshape input codeword to decode to N_cw rows
                decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
                % scrambler output stream is broken into blocks of L_cwd bits
                L_cwd = L_ldpc*MCS.CR;
                
                % preallocate output of the decoder
                decoderOut_cw_tmp = zeros(N_cw, L_cwd); % all cw
                
                % decode each code word
                for i_cw = 1:N_cw
                    decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                    
                    decoderOut_cw = hDec.step(decoderIn_cw_temp.');
                    decoderOut_cw_tmp(i_cw, :) = decoderOut_cw.';
                end
                
            elseif MCS.CR == 7/8
                % in this case, 48 bits are punctured from the parity bits of
                % the code rate 13/16
                % scrambler output stream is broken into blocks of L_cr78 bits
                L_cr78 = 546;
                N_remove_parity = 48; % number of parity bits to remove
                % preallocate output of the decoder
                decoderOut_cw_tmp = zeros(N_cw, L_cr78); % all cw
                
                % reshape input codeword to decode to N_cw rows
                decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
                
                % decode each code word using LDPC matrix for code rate 13/16
                for i_cw = 1:N_cw
                    decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                    % pick data part
                    cw_data_part = decoderIn_cw_temp(1, 1:L_cr78);
                    % pick parity part and add zeros (add first 48 parity bits
                    % to each codeword parity (logical 0 = +1))
                    cw_parity_part = [ones(1, N_remove_parity), decoderIn_cw_temp(1, L_cr78+1:end)];
                    
                    decoderIn_cw_temp2 = [cw_data_part, cw_parity_part]; % decoder input
                    
                    decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                    decoderOut_cw_tmp(i_cw, :) = decoderOut_cw.';
                end
                
            end
        case 2
            L_z = L_ldpc/(2*MCS.Repetition); % Length of block in bits
            L_zeros = L_z; % number of zeros
            
            % preallocate output of the decoder
            decoderOut_cw_tmp = zeros(N_cw, L_z); % all cw
            
            
            % reshape input codeword to decode to N_cw rows
            decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
            
            for i_cw = 1:N_cw
                decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                useRepeat = 1;
                % extract repetition (L_z bits)
                % for now, it is not employed
                decoderIn_cw_repeated_part_scrambled = decoderIn_cw_temp(1, L_z+1:L_z+L_zeros);
                PNseq_for_repetition = AD_ScramblingSeqGen(length(decoderIn_cw_repeated_part_scrambled), [1 1 1 1 1 1 1]);
                decoderIn_cw_repeated_part = ((decoderIn_cw_repeated_part_scrambled.').*(-2*PNseq_for_repetition+1)).';
                % replace repeated part in decoderIn_cw_temp by zero symbols (log. 0 = +1)
                decoderIn_cw_temp(1, L_z+1:L_z+L_zeros) = +10;
                
                if useRepeat == 1
                    decoderIn_cw_temp(1,1:L_z) = decoderIn_cw_temp(1,1:L_z)+decoderIn_cw_repeated_part;
                end
                decoderOut_cw_with_decoded_zeros = hDec.step(decoderIn_cw_temp.');
                decoderOut_cw_tmp(i_cw, :) = decoderOut_cw_with_decoded_zeros(1:L_z).';
            end
    end
    
    % build single decoder output stream
    decoderOut = reshape(decoderOut_cw_tmp.',[],1).';
    
elseif strcmp(MCS.phy_type,'OFDM') == 1 % OFDM -------------------------------------------------------------
    decoderIn_data_blocks = reshape(decoderIn,[],N_cw).';
    decoderOut_all_cw = zeros(N_cw, num_data_bits); % preallocation
    for i_cw = 1:N_cw
        decoderIn_cw_temp = decoderIn_data_blocks(i_cw, :);
        
        decoderOut_cw_temp = hDec.step(decoderIn_cw_temp.');
        decoderOut_cw_temp = decoderOut_cw_temp.';
        % save input word and parity to matrix (each row = 1 codeword)
        decoderOut_all_cw(i_cw, :) = double(decoderOut_cw_temp);
    end
    
    % build single decoder output stream
    decoderOut = reshape(decoderOut_all_cw.',[],1).';
    decoderOut_padded = decoderOut;
    
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC -------------------------------------------------------------
    % output from demodulator
    deintrlIn = double(demodulatorOut < 0); % slicer (hard)
    % de-interleaving
    deintrl_rows = 7;
    deintrl_cols = 8;
    N_deintrl = numel(deintrlIn)/(deintrl_rows*deintrl_cols);
    deintrlOut = cell(N_deintrl, 1);
    for i_deintrl = 1:N_deintrl
        deintrlIn_tmp = deintrlIn((i_deintrl-1)*56+1:i_deintrl*56,:);
        deintrlIn_tmp2 = reshape(deintrlIn_tmp, deintrl_rows, []).';
        deintrlOut{i_deintrl, 1} = deintrlIn_tmp2(:);
    end
    deintrlOut = vertcat(deintrlOut{:});
    % separate encoded data and scrambled pad bits
    blckIn = deintrlOut(1:end-N_blk_pad);
    % inner decoding (Block)
    if blck_n ~= blck_k
        blckOut = decode(blckIn, blck_n, blck_k, 'linear/binary', genmat, wifi_params.coding.Blck.SyndTable); % decoding
    else
        blckOut = blckIn;
    end
    % outer decoding (Reed Solomon)
    % broke data into blocks of length 208*8 bits
    rs_decoderIn_cw = cell(1, N_cw);
    rs_decoderOut_cw = cell(1, N_cw);
    
    decoderIn = blckOut;
    
    for i_cw = 1:N_cw
        if length(decoderIn) >= (208+16)*8
            rs_decoderIn_tmp = decoderIn(1:(208+16)*8);
            decoderIn(1:(208+16)*8) = [];
            rs_decoderIn_tmp_dec = bi2de(reshape(rs_decoderIn_tmp,8,[]).');
            rs_decoderIn_cw{1, i_cw} = rs_decoderIn_tmp_dec;
            % RS block decoding - all codewords except the last codeword
            rs_decoderOut_cw{1, i_cw} = step(hDec_RS{1}, rs_decoderIn_tmp_dec);
        else
            rs_decoderIn_tmp = decoderIn(1:end);
            decoderIn(1:end) = [];
            rs_decoderIn_tmp_dec = bi2de(reshape(rs_decoderIn_tmp,8,[]).');
            rs_decoderIn_cw{1, i_cw} = rs_decoderIn_tmp_dec;
            % RS block decoding - the last codeword
            rs_decoderOut_cw{1, i_cw} = step(hDec_RS{2}, rs_decoderIn_tmp_dec);
        end
        rs_decoderOut_dec = vertcat(rs_decoderOut_cw{:});
        rs_decoderOut_bin = de2bi(rs_decoderOut_dec,'left-msb').';
        rs_decoderOut = rs_decoderOut_bin(:);
        
    end
    
else
    
end
%% DESCRAMBLING
if strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    % descrambling sequence defined in line 100
    descramblerIn = decoderOut;
    if wifi_params.general.useScrambling == 0
        descramblerOut = xor(descramblerIn, zeros(size(scramblingSeq_data.')));
    else
        descramblerOut = xor(descramblerIn, scramblingSeq_data.');
    end
    
elseif strcmp(MCS.phy_type, 'OFDM') == 1 % OFDM --------------------------------------------------------
    descramblerIn = decoderOut;
    
    scr_seed = wifi_params.scrambling.scr_seed;
    scramblingSeq = AD_ScramblingSeqGen(length(descramblerIn)+504, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length(descramblerIn)).'; % scrambling seq for input data
    if wifi_params.general.useScrambling == 0
        descramblerOut = xor(descramblerIn, zeros(size(scramblingSeq_data))); % tx scrambling OFF
    else
        descramblerOut = xor(descramblerIn, scramblingSeq_data); % tx scrambling ON
    end
    
elseif strcmp(MCS.phy_type, 'LPSC') == 1 % LPSC --------------------------------------------------------
    scr_seed = wifi_params.scrambling.scr_seed;
    descramblingSeq = AD_ScramblingSeqGen(N_eb+N_blk_pad, scr_seed);
    descramblingSeq_data = descramblingSeq(1:length(rs_decoderOut)); % descrambling seq for input data
    
    if wifi_params.general.useScrambling == 0
        descramblerOut = rs_decoderOut; % tx scrambling OFF
    else
        descramblerOut = xor(rs_decoderOut, descramblingSeq_data); % tx scrambling ON
    end
    
else % Ctrl, OFDM -----------------------------------------------------------
    descramblerOut = decoderOut;
end
%% EXTRACT PSDU INFORMATION
% PSDU
if (strcmp(MCS.phy_type,'SC') == 1) || (strcmp(MCS.phy_type,'OFDM') == 1) % SC or OFDM ---------------------------------------------------------
    PSDU_rx_padded = descramblerOut;
    N_PSDU_padding_zeros = wifi_params.coding.N_Data_pad;
    
    PSDU_rx = PSDU_rx_padded(1,1:LENGTH*8);
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC ---------------------------------------------------------
    PSDU_rx_padded = descramblerOut.';
    PSDU_rx = descramblerOut.';
    
else
    % Ctrl -----------------------------------------------------------
    PSDU_rx_padded = decoderOut_padded;
    PSDU_rx = descramblerOut(1,5*8+1:end);
    Header_rx = descramblerOut(1, 1:5*8);
end
%% Save results

rxObj.output.PSDU = PSDU_rx.'; % user traffic data without overhead
rxObj.output.PSDU_padded = PSDU_rx_padded.'; % user traffic data without overhead padded with zeros

if strcmp(MCS.phy_type,'SC') == 1 % SC -----------------------------------------------------------------
    rxObj.output.coded_data = decoderOut.'; % decoder output - user traffic data + SERVICE zeros field with overhead
    rxObj.output.uncoded_data = (demodulatorOut(1:length(decoderIn)-N_blk_pad)<0).' + 0; % decoder input (raw data) (demodulated bits)
elseif strcmp(MCS.phy_type,'LPSC') == 1 % LPSC ---------------------------------------------------------
    rxObj.output.coded_data = rs_decoderOut.'; % decoder output - user traffic data + SERVICE zeros field with overhead
    rxObj.output.uncoded_data = deintrlIn.'; % decoder input (raw data) (demodulated bits)
elseif strcmp(MCS.phy_type,'OFDM') == 1
    rxObj.output.coded_data = decoderOut_padded.'; % decoder output - same as decoderOut in this case
    rxObj.output.uncoded_data = demodulatorOutHD_OFDM + 0; % decoder input (raw data) (demodulated bits)
else
    rxObj.output.coded_data = decoderOut_padded.'; % decoder output - user traffic data + SERVICE zeros field with overhead
    rxObj.output.uncoded_data = demodulatorOutHD_Ctrl.' + 0; % decoder input (raw data) (demodulated bits)
end
end