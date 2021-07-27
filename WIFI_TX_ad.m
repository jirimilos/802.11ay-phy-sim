function txObj = WIFI_TX_ad(wifi_params,txObj,mcs_num)
% WIFI transmitter function according to IEEE 802.11ad standard
% Author:   Jiri Milos, DREL FEEC BUT, 2019

%% DEFINE MODEL PARAMETERS
% Used MCS
MCS = wifi_params.MCS(mcs_num);
% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;
% Channel coding
switch MCS.coding_type
    case {'LDPC'}
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = wifi_params.coding.N_blks;
        % number of symbol padding bits
        N_blk_pad = wifi_params.coding.N_blk_pad;
        % LDPC encoder definition
        hEnc = comm.LDPCEncoder(wifi_params.coding.UsedLDPC_matrix);
    case {'RS+Blck', 'RS+SPC'}
        % number of RS codewords
        N_cw = wifi_params.coding.RS.N_cw;
        N_encsym_tot = wifi_params.coding.RS.N_encsym_tot;
        n_rs_encoders = wifi_params.coding.RS.numEncoders;
        hEnc_RS = cell(1, n_rs_encoders);
        % RS encoder(s) definition
        for ienc = 1:n_rs_encoders
            n = wifi_params.coding.RS.n_k(ienc, 1);
            k = wifi_params.coding.RS.n_k(ienc, 2);
            hEnc_RS{ienc} = comm.RSEncoder(...
                'BitInput', false,...
                'CodewordLength', n,...
                'MessageLength', k,...
                'PrimitivePolynomialSource', 'Property',...
                'PrimitivePolynomial', wifi_params.coding.RS.P_bin,...
                'GeneratorPolynomialSource','Property',...
                'GeneratorPolynomial', wifi_params.coding.RS.G{ienc});
        end
        
        % Block encoder definition
        blck_genmat = wifi_params.coding.Blck.G;
        blck_n = wifi_params.coding.Blck.n_k(1);
        blck_k = wifi_params.coding.Blck.n_k(2);
        N_eb = wifi_params.coding.Blck.N_eb;
        N_blk_pad = wifi_params.coding.Blck.N_blk_pad;
end

% Modulation mapping
hMod = wifi_params.modulation(log2(MCS.M)).hMod;
modNorm = wifi_params.modulation(log2(MCS.M)).hMod_factor;
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
if wifi_params.general.sendAllZeros == true
    PSDU_tx = zeros(LENGTH*8,1);
elseif wifi_params.general.sendAllOnes == true
    PSDU_tx = ones(LENGTH*8,1);
else
    PSDU_tx = randi([0 1],LENGTH*8,1);
end

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
else
    error('Define the rest IEEE 802.11ad PHY types.');
    
end
% Ctrl: Scrambling defined in Sec. 20.4.3.3.2, page 2454
% SC: both should be scrambled, padding zeros with continuous scrambling word, see page 2474
%% SCRAMBLING
if (strcmp(MCS.phy_type, 'Ctrl') == 1)
    scramblerOut = PSDU_tx_padded;
    scramblingSeq_data = [];
elseif (strcmp(MCS.phy_type, 'SC') == 1)
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
elseif (strcmp(MCS.phy_type, 'LPSC') == 1)
    scr_seed = wifi_params.scrambling.scr_seed;
    scramblingSeq = AD_ScramblingSeqGen(N_eb+N_blk_pad, scr_seed);
    scramblingSeq_data = scramblingSeq(1:length(PSDU_tx_padded)); % scrambling seq for input data
    scramblingSeq_block_pad = scramblingSeq(length(PSDU_tx_padded)+1:end); % scrambling sequence for blocks pad bits, see ch. 20.7.2.3.3.1/5, page 2482
    if wifi_params.general.useScrambling == 0
        scramblerOut = PSDU_tx_padded; % tx scrambling OFF
    else
        scramblerOut = xor(PSDU_tx_padded, scramblingSeq_data); % tx scrambling ON
    end
end

% scramblerOut = PSDU_tx_padded; % vyrazeni scrambleru
%% LDPC or RS+Blck CHANNEL CODING
encoderIn = scramblerOut;

if strcmp(MCS.phy_type,'Ctrl') == 1
    % preallocation of output auxiliary matrix
    encoderOut_all_cw = cell(1, N_cw);
    % scrambler output stream is broken into blocks of L_cwd bits
    L_cwd = 168;
    % sequences of length L_cwd, each row = single input data word
    encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
    encoderIn_data_blocks = [encoderIn_data_blocks, zeros(N_cw, (L_ldpc*3/4)-L_cwd)];
    
    % encode each data word
    for i_cw = 1:N_cw
        encoderIn_cw_temp = encoderIn_data_blocks(i_cw, :);
        
        encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
        encoderOut_cw_temp = encoderOut_cw_temp.';
        
        inputBits = encoderIn_data_blocks(i_cw, 1:168);
        parityBits = encoderOut_cw_temp(1, (L_ldpc*3/4)+1:end);
        % save input word and parity to matrix (each row = 1 codeword)
        encoderOut_all_cw{1, i_cw} = [inputBits(1:end-wifi_params.coding.N_Data_pad(i_cw)), parityBits]; % all zeros removed, including zeros defined by .coding.N_Data_pad (to have CR <= 1/2), also number of encoded bits is divisible by 32
        %         heatmap(encoderOut_all_cw) % show discrete values in 2D
    end
    % reshape to single codeword stream [1x?]
    encoderOut = [encoderOut_all_cw{:}];
    
    encoderOut_single_row = encoderOut;
elseif strcmp(MCS.phy_type,'SC') == 1
    
    % LDPC Encoding process is fully described in 20.6.3.2.3.3
    
    % a] is done in load_wifi_params.m and LDPC_params.m
    % b] procedure for converting the scrambled PSDU data to LDPC codeword
    switch MCS.Repetition
        case 1
            if MCS.CR ~= 7/8
                % preallocation of output auxiliary matrix
                encoderOut_all_cw = zeros(N_cw, L_ldpc);
                % scrambler output stream is broken into blocks of L_cwd bits
                L_cwd = L_ldpc*MCS.CR;
                % sequences of length L_cwd, each row = single input data word
                encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
                % encode each data word
                for i_cw = 1:N_cw
                    encoderIn_cw_temp = encoderIn_data_blocks(i_cw, :);
                    
                    encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                    encoderOut_cw_temp = encoderOut_cw_temp.';
                    % save input word and parity to matrix (each row = 1 codeword)
                    encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp;
                end
                
            elseif MCS.CR == 7/8
                % in this case, 48 bits are punctured from the parity bits of
                % the code rate 13/16
                % scrambler output stream is broken into blocks of L_cr78 bits
                L_cr78 = 546;
                N_remove_parity = 48; % number of parity bits to remove
                % preallocation of output auxiliary matrix
                encoderOut_all_cw = zeros(N_cw, L_ldpc+N_remove_parity);
                % sequences of length L_cr78, each row = single input data word
                encoderIn_data_blocks = reshape(encoderIn, L_cr78, []).';
                % encode each data word using LDPC matrix for code rate 13/16
                for i_cw = 1:N_cw
                    encoderIn_cw_temp = encoderIn_data_blocks(i_cw, :);
                    
                    encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                    encoderOut_cw_temp = encoderOut_cw_temp.';
                    % save input word and parity to matrix (each row = 1 codeword)
                    encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp;
                end
                % remove the first 48 parity bits from each codeword
                encoderOut_all_cw(:,L_cr78+1:L_cr78+N_remove_parity) = [];
                
            end
        case 2
            % preallocation of output auxiliary matrix
            encoderOut_all_cw = zeros(N_cw, L_ldpc);
            
            L_z = L_ldpc/(2*MCS.Repetition); % Length of block in bits
            L_zeros = L_z; % number of zeros
            % sequence length 2*L_z
            encoderIn_data_blocks = reshape(encoderIn, L_z, []).';
            encoderIn_data_blocks = [encoderIn_data_blocks, zeros(length(encoderIn)/L_z, L_zeros)]; % concatenate with zeros
            
            for i_cw = 1:N_cw
                encoderIn_cw_temp = encoderIn_data_blocks(i_cw, :);
                
                encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                encoderOut_cw_temp = encoderOut_cw_temp.';
                % repetiton and scrambling
                encoderIn_cw_for_repetition = encoderIn_cw_temp(1:L_z);
                PNseq_for_repetition = AD_ScramblingSeqGen(length(encoderIn_cw_for_repetition), [1 1 1 1 1 1 1]);
                encoderIn_cw_repeated_scrambled = xor(encoderIn_cw_for_repetition, PNseq_for_repetition.');
                % replacement of zeros bits (from L_z+1 to 336) by repeated and scrambled input sequence
                encoderOut_cw_temp2 = encoderOut_cw_temp;
                encoderOut_cw_temp2(1,L_z+1:336) = encoderIn_cw_repeated_scrambled;
                encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp2;
            end
            
        otherwise
            error('Wrong repetition factor, see WIFI_TX_ad.m and mcs_definition.m');
    end
    % c] concatenate to single codeword stream [1x?]
    encoderOut_single_row = reshape(encoderOut_all_cw.', 1, []);
    
    % d] divide to symbol blocks and add symbol pad bits
    % e] concatenation of coded bit stream and N_blk_pad zeros
    block_pad_zeros = zeros(N_blk_pad, 1);
    % scrambling of block pad bits with continuous scrambling sequence from
    % data scrambling
    if wifi_params.general.useScrambling == 0
        scrambled_block_pad_zeros = xor(block_pad_zeros, zeros(size(scramblingSeq_block_pad))); % tx block pad zeros scrambling OFF
    else
        scrambled_block_pad_zeros = xor(block_pad_zeros, scramblingSeq_block_pad); % tx block pad zeros scrambling ON
    end
    % concatenation
    encoderOut = [encoderOut_single_row, scrambled_block_pad_zeros.']; % for raw BER computing
    % symbol blocks (each row = 1 symbol block [bits])
    % symbol_blocks_bits = reshape(encoderOut.',[],N_blks).';
elseif strcmp(MCS.phy_type,'OFDM') == 1
    % preallocation of output auxiliary matrix
    encoderOut_all_cw = zeros(N_cw, L_ldpc);
    % scrambler output stream is broken into blocks of L_cwd bits
    L_cwd = L_ldpc*MCS.CR;
    % sequences of length L_cwd, each row = single input data word
    encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
    % encode each data word
    for i_cw = 1:N_cw
        encoderIn_cw_temp = encoderIn_data_blocks(i_cw, :);
        
        encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
        encoderOut_cw_temp = encoderOut_cw_temp.';
        % save input word and parity to matrix (each row = 1 codeword)
        encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp;
    end
    % concatenate to single codeword stream [1x?]
    encoderOut = reshape(encoderOut_all_cw.', 1, []);
    encoderOut_single_row = encoderOut;
    
elseif strcmp(MCS.phy_type,'LPSC') == 1
    % outer coding (Reed-Solomon)
    encoderIn_tmp = encoderIn;
    % broke data into blocks of length 208*8 bits
    rs_encoderIn_cw = cell(1, N_cw);
    rs_encoderOut_cw = cell(1, N_cw);
    for i_cw = 1:N_cw
        if length(encoderIn) >= 208*8
            rs_encoderIn_tmp = encoderIn(1:208*8);
            encoderIn(1:208*8) = [];
            rs_encoderIn_tmp_dec = bi2de(reshape(rs_encoderIn_tmp,8,[]).');
            rs_encoderIn_cw{1, i_cw} = rs_encoderIn_tmp_dec;
            % RS block encoding - all codewords except the last codeword
            rs_encoderOut_cw{1, i_cw} = step(hEnc_RS{1}, rs_encoderIn_tmp_dec);
        else
            rs_encoderIn_tmp = encoderIn(1:end);
            encoderIn(1:end) = [];
            rs_encoderIn_tmp_dec = bi2de(reshape(rs_encoderIn_tmp,8,[]).');
            rs_encoderIn_cw{1, i_cw} = rs_encoderIn_tmp_dec;
            % RS block encoding - the last codeword
            rs_encoderOut_cw{1, i_cw} = step(hEnc_RS{2}, rs_encoderIn_tmp_dec);
        end
    end
    rs_encoderOut_dec = vertcat(rs_encoderOut_cw{:});
    rs_encoderOut_bin = de2bi(rs_encoderOut_dec,'left-msb').';
    rs_encoderOut = rs_encoderOut_bin(:);
    % inner coding (Block)
    if blck_n ~= blck_k
        blckOut = encode(rs_encoderOut,blck_n,blck_k,'linear/binary',blck_genmat);
    else
        blckOut = rs_encoderOut;
    end
    % padding with zeros
    blckOut_pad = [blckOut; xor(scramblingSeq_block_pad(1:wifi_params.coding.Blck.N_blk_pad,:), zeros(wifi_params.coding.Blck.N_blk_pad,1))]; % should be scrambled, see 20.7.2.3.3.1/5
    % interleaving
    intrl_rows = 7;
    intrl_cols = 8;
    N_intrl = length(blckOut_pad)/(intrl_rows*intrl_cols);
    intrlIn = reshape(blckOut_pad, 8, []).';
    intrlOut = cell(N_intrl, 1);
    for i_intrl = 1:N_intrl
        intrlIn_tmp = intrlIn((i_intrl-1)*7+1:i_intrl*7,:);
        intrlOut{i_intrl, 1} = intrlIn_tmp(:);
    end
    intrlOut = vertcat(intrlOut{:});
    
    encoderOut = intrlOut.';
    encoderOut_single_row = encoderOut;
    
    encoderIn = encoderIn_tmp.';
else
    error('Unknown channel coder.');
end


%% MODULATION MAPPING
switch MCS.phy_type
    case 'Ctrl'
        % nondifferential stream, Sec. 20.4.3.3.4
        s_k = 2*encoderOut-1;
        
        d_k = zeros(size(s_k));
        for i_k = 1:length(s_k)
            if i_k == 1
                d_km1 = 1;
            else
                d_km1 = d_k(i_k-1);
            end
            d_k(1, i_k) = s_k(1, i_k)*d_km1;
        end
        % spreading, Sec. 20.4.3.3.5
        n = 0:(length(d_k)*32)-1;
        n_mod32 = mod(n, 32);
        N_Spr_Blks = length(n)/32;
        Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
        Ga32_spr = Ga32(n_mod32+1);
        i_n_flr_div32 = floor(n/32);
        d_k_spr = d_k(i_n_flr_div32+1);
        modulatedDataSymbol = Ga32_spr.*d_k_spr;
        modulatedDataSymbol_rotated = modulatedDataSymbol.*exp(1j*pi*n/2);
        
    case {'SC', 'LPSC'}
        % common modulation
        modulatedDataSymbol = step(hMod,encoderOut.')*modNorm;
        modulated_STF_symbol_index_k = (0:length(modulatedDataSymbol)-1).';
        % add pi/2 rotation
        if wifi_params.general.dataSymbolRotation == 0
            modulatedDataSymbol_rotated = modulatedDataSymbol;
        else
            modulatedDataSymbol_rotated = modulatedDataSymbol.*exp(1j*pi*modulated_STF_symbol_index_k/2);
        end
        
    case 'OFDM'
        N_mod_groups = numel(encoderOut)/MCS.N_cbps;
        % broke the input stream into groups of N_cbps bits
        modulatorInputGroups = reshape(encoderOut,[],N_mod_groups).';
        modulatedDataSymbolsGroups = zeros(N_mod_groups, wifi_params.mapping.n_data);
        
        % Spreaded QPSK (SQPSK)
        if (strcmp(MCS.MCS, 'MCS13') == 1) || (strcmp(MCS.MCS, 'MCS14') == 1) % SQPSK is used
            % Static Tone Pairing (STP) only
            for i_gr = 1:N_mod_groups
                modIputGrpTmp =  modulatorInputGroups(i_gr, :);
                % first half of the symbol modulation
                c_tmp = reshape(modIputGrpTmp,2,[]);
                dk_FH = (1/sqrt(2))*((2*c_tmp(1,:)-1)+1j*(2*c_tmp(2,:)-1));
                % second half of the symbol modulation
                dc_SH = conj(dk_FH);
                % concatenate first- and second half of modulated symbols -> SQPSK
                d_tmp = [dk_FH, dc_SH];
                %
                modulatedDataSymbolsGroups(i_gr,:) = d_tmp;
            end
            modulatedDataSymbol = reshape(modulatedDataSymbolsGroups.',[],1); % reshape into single stream
        else % QPSK, 16QAM and 64QAM
            for i_gr = 1:N_mod_groups % per mod. groups
                c_q =  modulatorInputGroups(i_gr, :);
                k = 168;
                Nsd = 336;
                % preallocation
                d_k = zeros(1,2*k);
                % choose groups of input bits, see Sec. 20.5.3.2.4.3 to 20.5.3.2.4.5
                if MCS.M == 4
                    c_qk = reshape(c_q, 2*MCS.N_bpsc, []);
                    for i_c = 1:k % per symbols 
                        d_out = DCM_mod(c_qk(:,i_c),'QPSK');
                        % Static tone pairing only (k = offset in frequency)
                        d_k(1, i_c) = d_out(1);
                        d_k(1, i_c+k) = d_out(2);
                    end
                    % stem(real(d_k)), drawnow                    
                elseif MCS.M == 16
                    % select two groups of 4 bits 
                    c_q1 = c_q(1:2*Nsd);
                    c_q2 = c_q(2*Nsd+1:length(c_q));
                    c_qk1 = reshape(c_q1, 4, []);
                    c_qk2 = reshape(c_q2, 4, []);
                    % concatenate it to single eigth of bits
                    c_qk = [c_qk1; c_qk2];
                    for i_c = 1:k % per symbols 
                        d_out = DCM_mod(c_qk(:,i_c),'16QAM');
                        % Static tone pairing only (k = offset in frequency)
                        d_k(1, i_c) = d_out(1);
                        d_k(1, i_c+k) = d_out(2);
                    end                    
                    % stem(real(d_k)), drawnow                    
                elseif MCS.M == 64
                    k64_vec = 0:Nsd-1;
                    m64_vec = (112*mod(k64_vec,3))+floor(k64_vec/3); % k mapped to m, as described in Sec. 20.5.3.2.4.5
                    % due to non DCM mapping in 64QAM modulation, the for loop
                    % has doubled number of iterations (for each subcarrier)
                    % interleaving between 3 LDPC codewords on a subcarrier basis
                    % reshape to groups of 6 bits
                    c_qk = reshape(c_q, 6, []);
                    for i_c = 1:Nsd % per symbols 
                        d_out = DCM_mod(c_qk(:,m64_vec(i_c)+1),'64QAM');
                        d_k(1,i_c) = d_out;
%                         stem(real(d_k)), drawnow
                    end
                else
                    error('Wrong DCM modulation type');
                end
                
                modulatedDataSymbolsGroups(i_gr,:) = d_k;    
%                 if i_gr == 1 % save genie information
%                     disp('temp DCM results TX saved, line 421');                    
%                     save('DCM.mat', 'c_q', 'c_qk','d_k');
%                 end
            end
            modulatedDataSymbol = reshape(modulatedDataSymbolsGroups.',[],1); % reshape into single stream
        end
    otherwise
        error('Unknown physical layer type.');
end

% konstelace(modulatedDataSymbol, modulatedDataSymbol_rotated)

%% SYMBOL BLOCKING AND GUARD INSERTION
if strcmp(MCS.phy_type,'Ctrl') == 1
    % NOP
    outputBlocks = modulatedDataSymbol_rotated;
    
elseif strcmp(MCS.phy_type,'SC') == 1
    % blocks division
    modulatedDataSymbol_blocks = reshape(modulatedDataSymbol_rotated,[],N_blks).';
    % guard insertion
    Ga64 = wifi_params.spreading.Golay_Seq.Ga_64;
    k_Ga64 = 0:n_gi-1;
    Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
    
    Ga64_rotated_repeated = repmat(Ga64_rotated, N_blks, 1);
    
    % concatenate guard and blocks
    modulatedDataSymbolsAndGuard = [Ga64_rotated_repeated, modulatedDataSymbol_blocks];
    % single stream of guard and blocks
    outputBlocks = reshape(modulatedDataSymbolsAndGuard.',[],1).';
    % add the last guard on the end of the final block
    outputBlocks = [outputBlocks, Ga64_rotated];
    
    % NOTE:
    % warning('Only data blocks are considered yet (see WIFI_TX_ad.m)');
elseif strcmp(MCS.phy_type,'OFDM') == 1
    % OFDM symbol creation
    % preallocation
    xOFDM = zeros(n_fft, wifi_params.coding.N_sym);
    % map pilot symbols, TBD
    mapPilots = repmat(wifi_params.mapping.map_pilots, 1, wifi_params.coding.N_sym);
    modulatedPilotSymbols = repmat(wifi_params.mapping.pilot_seq, 1, wifi_params.coding.N_sym); % see page 2465!
    xOFDM(mapPilots) = modulatedPilotSymbols;
    
    % map data symbols
    mapData = repmat(wifi_params.mapping.map_data, 1, wifi_params.coding.N_sym);
    xOFDM(mapData) = modulatedDataSymbol.';
    % perform IFFT
    N_tones = wifi_params.mapping.n_tot.nonzero; % number of active subcarriers
%     sqrt(n_fft)*sqrt(n_fft/wifi_params.mapping.n_tot.all)
%     xOFDM_ifft = (n_fft/sqrt(wifi_params.mapping.n_tot.all))*ifft(xOFDM, n_fft);
%     xOFDM_ifft = (n_fft/sqrt(N_tones))*ifft(xOFDM, n_fft); % see Sec. 20.5.3.2.6
%     xOFDM_ifft = sqrt(N_tones)*ifft(xOFDM, n_fft); % see Sec. 20.5.3.2.6 % unit power powA = (1/(2*512+1))*sum(abs(a.^2)) ??
    xOFDM_ifft = (n_fft/N_tones)*ifft(xOFDM, n_fft); % see Sec. 20.5.3.2.6
    % get cyclic prefix
    cyclic_prefix = xOFDM_ifft(wifi_params.cprefix.indices,:);
    % add cyclic prefix
    xOFDM_cp_ifft = [cyclic_prefix; xOFDM_ifft];
    % OFDM signal in single stream
    outputBlocks = reshape(xOFDM_cp_ifft,1,[]);
    
elseif strcmp(MCS.phy_type,'LPSC') == 1
    % see 20.7.2.3.5
    % guard insertion
    Ga64 = wifi_params.spreading.Golay_Seq.Ga_64;
    k_Ga64 = 0:n_gi-1;
    Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
    G8rotated = Ga64_rotated(end-7:end); % copy of the last 8 samples of Ga64
    % reshape modulated and rotated symbols to matrix with 56 columns
    N_d56_blocks = numel(modulatedDataSymbol_rotated)/56;
    outputBlocks = zeros(wifi_params.coding.Blck.N_blks*8,64);
    map_G64 = false(size(outputBlocks));
    map_G64((((1:wifi_params.coding.Blck.N_blks)-1)*8)+1, :) = true;
    map_d56andG8 = ~map_G64;
    d56_blocks = reshape(modulatedDataSymbol_rotated,56,[]).';
    % add G8 at the end of
    d56_blocks_G8 = [d56_blocks, repmat(G8rotated,N_d56_blocks,1)];
    % create blocks of 512 symbols
    outputBlocks(map_G64) = repmat(Ga64_rotated, wifi_params.coding.Blck.N_blks,1);
    outputBlocks(map_d56andG8) = d56_blocks_G8;
    % add the last Ga64 sequence
    outputBlocks = [outputBlocks; Ga64_rotated];
    outputBlocks = reshape(outputBlocks.',1,[]);
else
    error('Unknown physical layer type.');
end

%% Process other frame fields and concatenate to FRAME
% auxiliary Golay sequences
Ga128 = wifi_params.spreading.Golay_Seq.Ga_128;
Gb128 = wifi_params.spreading.Golay_Seq.Gb_128;
Gu512 = [-Gb128, -Ga128, +Gb128, -Ga128];
Gv512 = [-Gb128, +Ga128, -Gb128, -Ga128];

% Gu512 = [Gb128, Ga128, Gb128, Ga128]; % test purpose
% Gv512 = [Gb128, Ga128, Gb128, Ga128]; % test purpose

Gv128 = -Gb128;
% =========================================================================
switch MCS.phy_type
    case 'Ctrl' % Control -----------------------------------------------------------------
        % STF {48x +Ga128, 1x -Ga128}:
        STF = [repmat(Gb128, 1, 48), -Gb128, -Ga128];
        modulated_STF_symbol_index_k = (0:length(STF)-1);
        % add pi/2 rotation
        modulatedSTFSymbol_rotated = STF.*exp(1j*pi*modulated_STF_symbol_index_k/2);
        % CEF {Gu512, Gv512, Gv128}:
        CEF = [Gu512, Gv512, Gv128];
        modulated_CEF_symbol_index_k = (0:length(CEF)-1);
        % add pi/2 rotation
        modulatedCEFSymbol_rotated = CEF.*exp(1j*pi*modulated_CEF_symbol_index_k/2);
        % Beamforming training field (BTF) - in future
        BTF = [];
        modulatedBTFSymbol_rotated = BTF; % Header.*exp(1j*pi*modulated_Header_symbol_index_k/2);
        
        % x_s = zeros(1, sum(wifi_params.framing.ChipLengths));
        x_s(wifi_params.framing.ChipMap{1}) = 1*modulatedSTFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{2}) = 1*modulatedCEFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{3}) = outputBlocks;
        x_s(wifi_params.framing.ChipMap{4}) = 0*modulatedBTFSymbol_rotated;
        % aux
        modulatedHeaderSymbol_rotated1 = [];
        modulatedHeaderSymbol_rotated2 = [];
        
    case {'SC','LPSC'}% SC or LPSC ----------------------------------------------------------------------
        % these waveforms must be checked whether they are corresponding to the
        % standard, see on page 2446 and 2447
        % STF {16x +Ga128, 1x -Ga128}:
        STF = [repmat(Ga128, 1, 16), -Ga128];
        modulated_STF_symbol_index_k = (0:length(STF)-1);
        % add pi/2 rotation
        modulatedSTFSymbol_rotated = STF.*exp(1j*pi*modulated_STF_symbol_index_k/2);
        % CEF {Gu512, Gv512, Gv128}:
        CEF = [Gu512, Gv512, Gv128];
        modulated_CEF_symbol_index_k = (0:length(CEF)-1);
        % add pi/2 rotation
        modulatedCEFSymbol_rotated = CEF.*exp(1j*pi*modulated_CEF_symbol_index_k/2);
        % Header {64 bits}:
        if MCS.field_extended == 1
            switch MCS.field_base
                case 6
                    Base_Length1 = floor(N_blks*4/3)*42;
                    Base_Length2 = floor(floor(N_blks*56/39)*68.25);
                case 7
                    Base_Length1 = floor(floor(N_blks*4/3)*52.5);
                    Base_Length2 = floor(floor(N_blks*8/3)*68.25);
                case 8
                    Base_Length1 = floor(N_blks*4/3)*63;
                    Base_Length2 = floor(floor(N_blks*112/39)*68.25);
                case 9
                    Base_Length1 = floor(floor(N_blks*4/3)*68.25);
                    Base_Length2 = N_blks*210;
                case 10
                    Base_Length1 = floor(N_blks*8/3)*42;
                    Base_Length2 = N_blks*252;
                case 11
                    Base_Length1 = floor(floor(N_blks*8/3)*52.5);
                    Base_Length2 = N_blks*273;
                case 12
                    Base_Length1 = floor(N_blks*8/3)*63;
                    Base_Length2 = floor(floor(N_blks*56/13)*68.25);
                otherwise
                    error('Impossible combination of MCS fields, see mcs_definition.m');
            end
            
            % Length_field in header
            Length_field = Base_Length1-floor((Base_Length2-wifi_params.general.LENGTH)/4);
        else
            % Length_field in header
            Length_field = wifi_params.general.LENGTH;
        end
        
        Length_field_bin = de2bi(Length_field, 18, 'right-msb');
        % (provisionally), see: 20.6.3.1.4 Header encoding and modulation
        % LSB first
        Header = [...
            wifi_params.scrambling.scr_seed,...         % bits 0...6 = scrambilng seed (length 7)
            de2bi(MCS.field_base,5,'right-msb'),... % bits 7...11 = MCS according to Table 20-17 (length 5)
            Length_field_bin,...                        % bits 12...29 = Length (length 18)
            0,...                                       % bits 30...30 = Additional PPDU (length 1) % NOTE: set to zero for now, see Table 20-17
            0,...                                       % bits 31...31 = Packet Type (length 1) % NOTE: set to zero for now, see Table 20-17
            [0 0 0 0 0],...                             % bits 32...36 = Training Length (length 5) % NOTE: set to all zeros vector for now - NO BEAMFORMING TRN, see Table 20-17
            0,...                                       % bits 37...37 = Aggregation (length 1) % NOTE: set to zero for now, see Table 20-17
            0,...                                       % bits 38...38 = Beam Tracking Request (length 1) % NOTE: set to zero for now, see Table 20-17
            [1 1 1 1],...                               % bits 39...42 = Last RSSI (length 4) % NOTE: set to all ones vector for now, see Table 20-17 (Last RSSI >= -42 dBm)
            0,...                                       % bits 43...43 = Turnaround (length 1) % NOTE: set to zero for now, see Table 20-17 and Table 20-1
            MCS.field_extended,...                      % bits 44...44 = Extended SC MCS Indication (length 1) % see Table 20-17
            [0 0 0],...                                 % bits 45...47 = Reserved all zero bits (length 3); shall be ignored in receiver % see Table 20-17
            ];
        HeaderCRC = CRC16_FCS(wifi_params.framing.Header.hCRCgen, Header).'; % Header CRC-16 calculation
        HeaderWithCRC = [Header, HeaderCRC];
        
        % header encoding
        % 1/ scrambling input header bits excepting first seven bits
        scramblingSeqHeader = AD_ScramblingSeqGen(length(HeaderWithCRC)-length(wifi_params.scrambling.scr_seed), wifi_params.scrambling.scr_seed).';
        % scrambling and concatenation with scrambling seq. results to $\bar{d}_{1s}=[q_1, q_2, \dots, q_{64}]$
        scramblerHeaderOut = [HeaderWithCRC(1:length(wifi_params.scrambling.scr_seed)), xor(HeaderWithCRC(length(wifi_params.scrambling.scr_seed)+1:end), scramblingSeqHeader)];
        % 2/ LDPC encoding (CR = 3/4)
        encoderIn_Header = [scramblerHeaderOut, zeros(1, 504-length(HeaderWithCRC))];
        
        hEncHeader = comm.LDPCEncoder(wifi_params.framing.Header.coding.UsedLDPC_matrix);
        encoderOut_Header = hEncHeader.step(encoderIn_Header.');
        % codeword sequence 1
        % length(cs1) = 224
        cs1 = encoderOut_Header;
        cs1([length(HeaderWithCRC)+1:504, 665:672]) = []; % remove defined bits
        % codeword sequence 2
        % length(cs2) = 224
        cs2 = encoderOut_Header;
        cs2([length(HeaderWithCRC)+1:504, 657:664]) = []; % remove defined bits
        
        scramblingSeq_cs2 = AD_ScramblingSeqGen(length(cs2), ones(1, 7)).';
        cs2_Scrambled = single(xor(cs2.', scramblingSeq_cs2));
        % concatenate cs1 and cs2_scrambled to form a 448 bits header codeword
        cs12 = [cs1.', cs2_Scrambled];
        
        % pi/2-BPSK mapping of header block
        modulated_Header_symbol_index_k = (0:length(cs12)-1);
        modulatedHeaderSymbol_rotated1 = (2*cs12-1).*exp(1j*pi*modulated_Header_symbol_index_k/2);
        modulatedHeaderSymbol_rotated_withGI1 = [Ga64_rotated, modulatedHeaderSymbol_rotated1]; % guard insertion
        % repeat HeaderSymbol block and multiply by -1
        modulatedHeaderSymbol_rotated2 = -1*modulatedHeaderSymbol_rotated1;
        modulatedHeaderSymbol_rotated_withGI2 = [Ga64_rotated, modulatedHeaderSymbol_rotated2];
        % concatenate modulated HeaderSymbol blocks with GI
        modulatedHeaderSymbol_rotated_withGI12 = [modulatedHeaderSymbol_rotated_withGI1, modulatedHeaderSymbol_rotated_withGI2];
        
        % Beamforming training field (BTF) - in future
        BTF = [];
        modulatedBTFSymbol_rotated = BTF; % Header.*exp(1j*pi*modulated_Header_symbol_index_k/2);
        
        % x_s = zeros(1, sum(wifi_params.framing.ChipLengths));
        x_s(wifi_params.framing.ChipMap{1}) = 1*modulatedSTFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{2}) = 1*modulatedCEFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{3}) = 1*modulatedHeaderSymbol_rotated_withGI12;
        x_s(wifi_params.framing.ChipMap{4}) = outputBlocks;
        x_s(wifi_params.framing.ChipMap{5}) = 0*modulatedBTFSymbol_rotated;
        
    case 'OFDM' % OFDM --------------------------------------------------------------------
        % windowing / only in OFDM (No windowing in CEF and STF)
        % these waveforms must be checked whether they are corresponding to the
        % standard, see on page 2446 and 2447
        
        % STF {16x +Ga128, 1x -Ga128}:
        STF = [repmat(Ga128, 1, 16), -Ga128];
        modulated_STF_symbol_index_k = (0:length(STF)-1);
        % add pi/2 rotation
        modulatedSTFSymbol_rotated = STF.*exp(1j*pi*modulated_STF_symbol_index_k/2);
        
        % CEF {Gv512, Gu512, Gv128}:
        CEF = [Gv512, Gu512, Gv128];
        modulated_CEF_symbol_index_k = (0:length(CEF)-1);
        % add pi/2 rotation
        modulatedCEFSymbol_rotated = CEF.*exp(1j*pi*modulated_CEF_symbol_index_k/2);
        % Header {64 bits}:
        % (provisionally), see: 20.6.3.1.4 Header encoding and modulation
        Length_field_bin = de2bi(wifi_params.general.LENGTH, 18, 'right-msb');
        % LSB first
        Header = [...
            wifi_params.scrambling.scr_seed,...         % bits 0...6 = scrambling seed (length 7)
            de2bi(MCS.field_base,5,'right-msb'),...     % bits 7...11 = MCS according to Table 20-13 (length 5)
            Length_field_bin,...                        % bits 12...29 = Length (length 18)
            0,...                                       % bits 30...30 = Additional PPDU (length 1) % NOTE: set to zero for now, see Table 20-13
            0,...                                       % bits 31...31 = Packet Type (length 1) % NOTE: set to zero for now, see Table 20-13
            [0 0 0 0 0],...                             % bits 32...36 = Training Length (length 5) % NOTE: set to all zeros vector for now - NO BEAMFORMING TRN, see Table 20-13
            0,...                                       % bits 37...37 = Aggregation (length 1) % NOTE: set to zero for now, see Table 20-13
            0,...                                       % bits 38...38 = Beam Tracking Request (length 1) % NOTE: set to zero for now, see Table 20-13
            0,...                                       % bits 39...39 = Tone Pairing Type (length 1) % NOTE: set to zero = STP; set to one = DTP, see Table 20-13
            0,...                                       % bits 40...40 = DTP Indicator (length 1) % NOTE: set to zero for now (reserved), see Table 20-13
            [0 0 0 0],...                               % bits 41...44 = Last RSSI (length 4) % NOTE: value 0, previous packet was not received, see Table 20-13
            0,...                                       % bits 45...45 = Turnaround (length 1) % NOTE: value 0, for now, see Table 20-13 and 20-1
            [0 0],...                                   % bits 46...47 = Reserved all zero bits (length 2); shall be ignored in receiver % see Table 20-13
            ];
        HeaderCRC = CRC16_FCS(wifi_params.framing.Header.hCRCgen, Header).'; % Header CRC-16 calculation
        HeaderWithCRC = [Header, HeaderCRC];
        % header encoding; header is encoded using 1 OFDM symbol
        % 1/ scrambling input header bits excepting first seven bits
        scramblingSeqHeader = AD_ScramblingSeqGen(length(HeaderWithCRC)-length(wifi_params.scrambling.scr_seed), wifi_params.scrambling.scr_seed).';
        % scrambling and concatenation with scrambling seq. results to $q=[q_1, q_2, \dots, q_{64}]$
        scramblerHeaderOut = [HeaderWithCRC(1:length(wifi_params.scrambling.scr_seed)), xor(HeaderWithCRC(length(wifi_params.scrambling.scr_seed)+1:end), scramblingSeqHeader)];
        % 2/ LDPC encoding (CR = 3/4)
        encoderIn_Header = [scramblerHeaderOut, zeros(1, 504-length(HeaderWithCRC))];
        
        hEncHeader = comm.LDPCEncoder(wifi_params.framing.Header.coding.UsedLDPC_matrix);
        encoderOut_Header = hEncHeader.step(encoderIn_Header.');
        encoderOut_Header_Parity = encoderOut_Header(end-168+1:end);
        % codeword sequence 1
        % length(cs1) = 64 scrambled header bits + 160 parity bits (index 9 to 168) = 224
        cs1 = [encoderOut_Header(1:64); encoderOut_Header_Parity(9:end)];
        % codeword sequence 2
        % length(cs2) = 64 + 84 + 76 = 224
        cs2 = [encoderOut_Header(1:64); encoderOut_Header_Parity(1:84); encoderOut_Header_Parity(93:end)];
        scramblingSeq_cs2 = AD_ScramblingSeqGen(length(cs2)+64+160, ones(1, 7)).';
        cs2_Scrambled = single(xor(cs2, scramblingSeq_cs2(1:length(cs2)).'));
        % codeword sequence 3
        % length(cs3) = 224
        cs3 = [encoderOut_Header(1:64); encoderOut_Header_Parity(1:160)];
        scramblingSeq_cs3 = scramblingSeq_cs2(length(cs2)+1:end);
        cs3_Scrambled = single(xor(cs3, scramblingSeq_cs3.'));
        % concatenate cs1 and cs2_scrambled to form a 448 bits header codeword
        cs123 = [cs1; cs2_Scrambled; cs3_Scrambled];
        
        % modulate Header input bits using QPSK (see Sec. 20.5.3.2.4.3)
        modulatorInputGroup_Header = cs123.';
        % broke into 4 bits groups
        cs123_tmp = reshape(modulatorInputGroup_Header,4,[]);
        x_cs123_1 = (1/sqrt(2))*((2*cs123_tmp(1,:)-1)+1j*(2*cs123_tmp(3,:)-1));
        x_cs123_2 = (1/sqrt(2))*((2*cs123_tmp(2,:)-1)+1j*(2*cs123_tmp(4,:)-1));
        x_cs123 = zeros(1,length(x_cs123_1)+length(x_cs123_2));
        x_cs123(1:2:end) = x_cs123_1;
        x_cs123(2:2:end) = x_cs123_2;
        % OFDM symbol creation
        % preallocation
        xOFDM_Header = zeros(n_fft, 1);
        
        % map pilot symbols
        mapPilots = wifi_params.mapping.map_pilots;
        modulatedPilotSymbols = wifi_params.mapping.pilot_seq;
        xOFDM_Header(mapPilots) = modulatedPilotSymbols;
        
        % map Header symbols
        mapData = wifi_params.mapping.map_data;
        xOFDM_Header(mapData) = x_cs123;
        
        % perform IFFT
        xOFDM_Header_ifft = (1/sqrt(N_tones))*ifft(xOFDM_Header, n_fft); % see Sec. 20.5.3.2.6
        % get cyclic prefix
        cyclic_prefix_H = xOFDM_Header_ifft(wifi_params.cprefix.indices,:);
        % add cyclic prefix
        xOFDM_Header_cp_ifft = [cyclic_prefix_H; xOFDM_Header_ifft];
        % Header OFDM signal in time domain
        modulatedHeaderSymbols = xOFDM_Header_cp_ifft.';
        
        % Beamforming training field (BTF) - in future
        BTF = [];
        modulatedBTFSymbol_rotated = BTF;
        
        % x_s = zeros(1, sum(wifi_params.framing.ChipLengths));
        x_s(wifi_params.framing.ChipMap{1}) = 1*modulatedSTFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{2}) = 1*modulatedCEFSymbol_rotated;
        x_s(wifi_params.framing.ChipMap{3}) = 1*modulatedHeaderSymbols;
        x_s(wifi_params.framing.ChipMap{4}) = outputBlocks;
        x_s(wifi_params.framing.ChipMap{5}) = 0*modulatedBTFSymbol_rotated;
        
    otherwise
        error('Undefined PHY type.')
end

%% Save outputs
txObj.output.PSDU = PSDU_tx; % user traffic data without overhead
txObj.output.PSDU_padded = PSDU_tx_padded; % user traffic data without overhead padded with zeros

txObj.output.coded_data = encoderIn; % encoder input
txObj.output.uncoded_data = encoderOut_single_row; % encoder output (raw data)

% txObj.genie.scramblingSeq.data = scramblingSeq_data; % scrambling sequence for data
% txObj.genie.scramblingSeq.block_pad = scramblingSeq_block_pad; % scrambling sequence for block pad bits

% save('tx_test.mat','PSDU_tx','PSDU_tx_plus_service','encoderIn','encoderOut');
txObj.output.Data = outputBlocks;
txObj.output.DataDerotated = modulatedDataSymbol.';
txObj.output.STF = STF;
txObj.output.CEF = CEF;
txObj.output.Header = Header;
txObj.output.HeaderCRC = HeaderCRC;
txObj.output.BTF = BTF;
txObj.output.x_s = x_s;

if strcmp(MCS.phy_type,'OFDM') == 0
    txObj.genie.modulatedDataSymbol = modulatedDataSymbol;
    txObj.genie.modulatedDataSymbol_rotated = modulatedDataSymbol_rotated;
    
    txObj.genie.x_header1 = modulatedHeaderSymbol_rotated1;
    txObj.genie.x_header2 = modulatedHeaderSymbol_rotated2;
else
    txObj.genie.modulatedDataSymbol = modulatedDataSymbol;
    txObj.genie.modulatedDataSymbol_rotated = [];
    
    txObj.genie.x_header = modulatedHeaderSymbols;
end

txObj.genie.x_STF = modulatedSTFSymbol_rotated;
txObj.genie.x_CEF = modulatedCEFSymbol_rotated;

txObj.genie.scramblingSeq_data = scramblingSeq_data;

end