function [decoderOut] = ayLDPCDecoding(decoderIn, wifi_params, MCS)

cbpb_index = 2; % normal length, TEST

SC_LDPC_softDec = true; % LDPC soft decoding indication
% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;

% normFaktor = [1 1/sqrt(2) 1 1/sqrt(10) -1 1/sqrt(42)];
normFaktor = [1 1 1 1/sqrt(10) -1 1/sqrt(42)];
% decoderIn = decoderIn*(1/wifi_params.modulation(log2(MCS.M)).hMod_factor);
decoderIn = decoderIn/normFaktor(log2(MCS.M));

if strcmp(MCS.phy_type,'Ctrl') == 1 % Control -----------------------------------------------------------------
    
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
    
    % number of codewords
    N_cw = wifi_params.coding.N_cw;
    % single codeword length
    L_ldpc = wifi_params.coding.L_ldpc;
    % number of symbol blocks
    N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
    % number of symbol padding bits
    N_blk_pad = (N_blks*MCS.N_cbpb(cbpb_index))-(N_cw*L_ldpc);
    % LDPC decoder definition
    if SC_LDPC_softDec == 1
        hDec = comm.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied','DecisionMethod','Soft decision'); % for soft decoding
    else
        hDec = comm.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied'); % for hard decoding
        % hDec = comm.gpu.LDPCDecoder(wifi_params.coding.UsedLDPC_matrix,'IterationTerminationCondition','Parity check satisfied');
    end
    
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
    if MCS.Repetition == 1 % no repetition =========================================================
        % divide according to codeword length:
        switch wifi_params.coding.L_ldpc
            case {672, 1344} % opt. a) --------------------------------
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
                
            case 624 % opt. b) ----------------------------------------
                % CR = 7/8
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
                    %                     cw_parity_part = [10*ones(1, N_remove_parity), decoderIn_cw_temp(1, L_cr78+1:end)];
                    cw_parity_part = [zeros(1, N_remove_parity), decoderIn_cw_temp(1, L_cr78+1:end)];
                    
                    decoderIn_cw_temp2 = [cw_data_part, cw_parity_part]; % decoder input
                    
                    decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                    decoderOut_cw_tmp(i_cw, :) = decoderOut_cw.';
                end
                
            case 1248 % opt. c) ---------------------------------------
                % CR = 7/8
                % in this case, 96 bits are punctured from the parity bits of
                % the code rate 13/16
                % scrambler output stream is broken into blocks of L_cr78 bits
                L_cr78 = 1092;
                N_remove_parity = 96; % number of parity bits to remove
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
                    cw_parity_part = [10*ones(1, N_remove_parity), decoderIn_cw_temp(1, L_cr78+1:end)];
                    
                    decoderIn_cw_temp2 = [cw_data_part, cw_parity_part]; % decoder input
                    
                    decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                    decoderOut_cw_tmp(i_cw, :) = decoderOut_cw.';
                end
                
            case {504, 1008} % opt. f), g), h) and i) -----------------
                if MCS.CR == 2/3 % 8PSK only, opt. f) and g) - - - - -
                    % reshape input codeword to decode to N_cw rows
                    decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
                    % scrambler output stream is broken into blocks of 336 bits
                    L_cwd = 336;
                    N_zero_bits = 168;
                    % preallocate output of the decoder
                    decoderOut_cw_tmp = zeros(N_cw, L_cwd); % all cw
                    
                    % decode each code word
                    for i_cw = 1:N_cw
                        decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                        decoderIn_cw_temp2 = [decoderIn_cw_temp(1:L_cwd), +10*ones(1, N_zero_bits), decoderIn_cw_temp(L_cwd+1:end)];
                        
                        decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                        decoderOut_cw_tmp(i_cw, :) = decoderOut_cw(1:L_cwd).';
                    end
                    
                elseif MCS.CR == 5/6 % 8PSK only, opt. h) and i) - - -
                    % reshape input codeword to decode to N_cw rows
                    decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
                    % scrambler output stream is broken into blocks of 420 bits
                    L_cwd = 420;
                    N_zero_bits = 168;
                    % preallocate output of the decoder
                    decoderOut_cw_tmp = zeros(N_cw, L_cwd); % all cw
                    
                    % decode each code word
                    for i_cw = 1:N_cw
                        decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                        decoderIn_cw_temp2 = [decoderIn_cw_temp(1:L_cwd), 10*ones(1, N_zero_bits), decoderIn_cw_temp(L_cwd+1:end)];
                        
                        decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                        decoderOut_cw_tmp(i_cw, :) = decoderOut_cw.';
                    end
                    
                else
                    error('AY: Wrong LDPC code rate');
                end
                
            case {468, 936} % opt j) and k) ---------------------------
                % reshape input codeword to decode to N_cw rows
                decoderIn_all_cw = reshape(decoderIn_extracted,[],N_cw).';
                % scrambler output stream is broken into blocks of 390 bits
                L_cwd = 390;
                N_zero_bits = 156;
                N_remove_parity = 48; % number of zero parity bits to add
                % preallocate output of the decoder
                decoderOut_cw_tmp = zeros(N_cw, L_cwd); % all cw
                
                % decode each code word
                for i_cw = 1:N_cw
                    decoderIn_cw_temp = decoderIn_all_cw(i_cw, :);
                    % pick data part
                    cw_data_part = [decoderIn_cw_temp(1, 1:L_cwd), 10*ones(1, N_zero_bits)];
                    % pick parity part and add zeros (add first 48 parity bits
                    % to each codeword parity (logical 0 = +1))
                    %                     cw_parity_part = [10*ones(1, N_remove_parity), decoderIn_cw_temp(1, L_cr78+1:end)];
                    cw_parity_part = [zeros(1, N_remove_parity), decoderIn_cw_temp(1, L_cwd+1:end)];
                    
                    decoderIn_cw_temp2 = [cw_data_part, cw_parity_part]; % decoder input
                    
                    decoderOut_cw = hDec.step(decoderIn_cw_temp2.');
                    decoderOut_cw_tmp(i_cw, :) = decoderOut_cw(1:L_cwd).';
                end
                
                
            otherwise
                error('AY: Wrong LDPC codeword length');
        end
        
    elseif MCS.Repetition == 2 % repetition ====================================================
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
            decoderIn_cw_temp(1, L_z+1:L_z+L_zeros) = +10; % strong zero (10 times)
            
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
    
end

if SC_LDPC_softDec == 1
    decoderOut = double(decoderOut<0).'; % for soft decoding
else
    decoderOut = decoderOut.'; % for hard decoding
end
end