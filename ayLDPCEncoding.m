function [encoderOut, encoderOut_single_row] = ayLDPCEncoding(encoderIn, wifi_params, MCS, scramblingSeq_block_pad)

% number of codewords
N_cw = wifi_params.coding.N_cw;
% single codeword length
L_ldpc = wifi_params.coding.L_ldpc;
% number of symbol padding bits
N_blk_pad = wifi_params.coding.N_blk_pad;
% LDPC encoder definition
hEnc = comm.LDPCEncoder(wifi_params.coding.UsedLDPC_matrix);

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
    % ad: LDPC Encoding process is fully described in 20.6.3.2.3.3
    % ay: LDPC Encoding process as described in IEEE 802.11-11/1712r00
    
    % a] is done in load_wifi_params.m and LDPC_params.m
    % b] procedure for converting the scrambled PSDU data to LDPC codeword
    if MCS.Repetition == 1 % no repetition ============================================
        % divide according to codeword length:
        switch wifi_params.coding.L_ldpc
            case {672, 1344} % opt. a) --------------------------------
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
                
            case 624 % opt. b) ----------------------------------------
                % CR = 7/8
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
                
            case 1248 % opt. c) ---------------------------------------
                % CR = 7/8
                % in this case, 96 bits are punctured from the parity bits of
                % the code rate 13/16
                % scrambler output stream is broken into blocks of L_cr78 bits
                L_cr78 = 1092;
                N_remove_parity = 96; % number of parity bits to remove
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
                
                
            case {504, 1008} % opt. f), g), h) and i) -----------------
                if MCS.CR == 2/3 % 8PSK only, opt. f) and g) - - - - -
                    % preallocation of output auxiliary matrix
                    encoderOut_all_cw = zeros(N_cw, L_ldpc);
                    % scrambler output stream is broken into blocks of 336 bits
                    L_cwd = 336;
                    N_zero_bits = 168;
                    % sequences of length L_cwd, each row = single input data word
                    encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
                    % encode each data word
                    for i_cw = 1:N_cw
                        encoderIn_cw_temp = [encoderIn_data_blocks(i_cw, :), false(1, N_zero_bits)]; % add 168 (=N_zero_bits) zero bits to create input word
                        
                        encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                        encoderOut_cw_temp = encoderOut_cw_temp.';
                        % save input word and parity to matrix (each row = 1 codeword)
                        % discard 168 above added zeros; keep 336 bits of input word and 168 bits of parity, only
                        encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp(1,[1:L_cwd, L_cwd+N_zero_bits+1:672]);
                    end
                    
                    
                elseif MCS.CR == 5/6 % 8PSK only, opt. h) and i) - - -
                    % preallocation of output auxiliary matrix
                    encoderOut_all_cw = zeros(N_cw, L_ldpc);
                    % scrambler output stream is broken into blocks of 420 bits
                    L_cwd = 420;
                    N_zero_bits = 168;
                    % sequences of length L_cwd, each row = single input data word
                    encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
                    % encode each data word
                    for i_cw = 1:N_cw
                        encoderIn_cw_temp = [encoderIn_data_blocks(i_cw, :), false(1, N_zero_bits)]; % add 168 (=N_zero_bits) zero bits to create input word
                        
                        encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                        encoderOut_cw_temp = encoderOut_cw_temp.';
                        % save input word and parity to matrix (each row = 1 codeword)
                        % discard 168 above added zeros; keep 336 bits of input word and 84 bits of parity, only
                        encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp(1,[1:L_cwd, L_cwd+N_zero_bits+1:672]);
                    end
                    
                else
                    error('AY: Wrong LDPC code rate');
                end
            case {468, 936} % opt j) and k) ---------------------------
                % preallocation of output auxiliary matrix
                encoderOut_all_cw = zeros(N_cw, L_ldpc);
                % scrambler output stream is broken into blocks of 390 bits
                L_cwd = 390;
                N_zero_bits = 156;
                N_remove_parity = 48; % number of parity bits to remove
                % sequences of length L_cwd, each row = single input data word
                encoderIn_data_blocks = reshape(encoderIn, L_cwd, []).';
                % encode each data word
                for i_cw = 1:N_cw
                    encoderIn_cw_temp = [encoderIn_data_blocks(i_cw, :), false(1, N_zero_bits)]; % add 156 (=N_zero_bits) zero bits to create input word
                    
                    encoderOut_cw_temp = hEnc.step(encoderIn_cw_temp.');
                    encoderOut_cw_temp = encoderOut_cw_temp.';
                    % save input word and parity to matrix (each row = 1 codeword)
                    % discard 168 above added zeros; keep 336 bits of input word and 84 bits of parity, only
                    encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp(1,[1:L_cwd, L_cwd+N_zero_bits+N_remove_parity+1:672]);
                end
                
            otherwise
                error('AY: Wrong LDPC codeword length');
        end
        
    elseif MCS.Repetition == 2 % repetition, MCS1 only ============================================
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
            encoderOut_cw_temp2(1,L_z+1:2*L_z) = encoderIn_cw_repeated_scrambled;
            encoderOut_all_cw(i_cw, :) = encoderOut_cw_temp2;
        end
        
    else
        error('AY: Wrong repetition factor, see WIFI_TX_ay.m and mcs_definition.m');
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
    
    
else
    error('Unknown channel coder.');
end

encoderOut = encoderOut.'; % column vector

end

