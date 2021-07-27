function [x_s, STF, CEF, Header, BTF] = ayBlockingAndGI(modulatedDataSymbol_rotated, wifi_params, MCS)

% auxiliary Golay sequences
Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
Gb32 = wifi_params.spreading.Golay_Seq.Gb_32;

Ga64 = wifi_params.spreading.Golay_Seq.Ga_64;
Gb64 = wifi_params.spreading.Golay_Seq.Gb_64;

Ga128 = wifi_params.spreading.Golay_Seq.Ga_128;
Gb128 = wifi_params.spreading.Golay_Seq.Gb_128;

Gu512 = [-Gb128, -Ga128, +Gb128, -Ga128];
Gv512 = [-Gb128, +Ga128, -Gb128, -Ga128];
Gv128 = -Gb128;

n_fft = wifi_params.mapping.n_fft; % all samples
n_gi = wifi_params.mapping.N_GI; % GI samples
N_SPB = wifi_params.mapping.N_SPB; % samples within a block
n_symb = N_SPB;

k_Ga32 = 0:32-1;
Ga32_rotated = Ga32.*exp(1j*pi*k_Ga32/2);
Gb32_rotated = Gb32.*exp(1j*pi*k_Ga32/2);

k_Ga64 = 0:64-1;
Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
Gb64_rotated = Gb64.*exp(1j*pi*k_Ga64/2);

k_Ga128 = 0:128-1;
Ga128_rotated = Ga128.*exp(1j*pi*k_Ga128/2);
Gb128_rotated = Gb128.*exp(1j*pi*k_Ga128/2);

nTxAnt = wifi_params.antConfig.nTxAnt;

if n_gi == 32
elseif n_gi == 64
else
end

if wifi_params.general.dataSymbolRotation == 1 % GI rotation
end

%% Processing of data blocks
% number of symbol blocks
N_blks = wifi_params.coding.N_blks;

if strcmp(MCS.phy_type,'Ctrl') == 1
    % NOP
    outputBlocks = modulatedDataSymbol_rotated;
    
elseif strcmp(MCS.phy_type,'SC') == 1
    
    switch wifi_params.antConfig.mode
        case {'SISO','RxD'} % =============================================
            %         modulatedDataSymbol_rotated = (1:length(modulatedDataSymbol_rotated)).';
            
            % blocks division
            modulatedDataSymbol_blocks = reshape(modulatedDataSymbol_rotated,[],N_blks).';
            % guard insertion (types: GI, CP or zeros)
            GIType = wifi_params.general.guardInterval.type;
            if strcmp(GIType,'GI') == 1
                Ga64_rotated_repeated = repmat(Ga64_rotated, N_blks, 1);
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [Ga64_rotated_repeated, modulatedDataSymbol_blocks];
                theLastGI = Ga64_rotated;
            elseif strcmp(GIType,'CP') == 1
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [modulatedDataSymbol_blocks(:,N_SPB-n_gi+1:N_SPB), modulatedDataSymbol_blocks];
                theLastGI = modulatedDataSymbol_blocks(N_blks,N_SPB-n_gi+1:N_SPB);
            else % zeros only
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [zeros(N_blks, n_gi), modulatedDataSymbol_blocks];
                theLastGI = zeros(1,n_gi);
            end
            % single stream of guard and blocks
            outputBlocks = reshape(modulatedDataSymbolsAndGuard.',[],1).';
            % add the last guard on the end of the final block
            outputBlocks = [outputBlocks, theLastGI];
            
            
        case 'TxD' % ======================================================
            % blocks division
            % Al Dhahir STBC Encoding: N_blks must be even!
            % Al Dhahir STBC Encoding + CP/GI:
            % in case of common GI (Ga64) the input of AlDhahir STBC Encoder
            % should be appended Data (448 samples)---Ga64 (64 samples)
            STBC_encoderInputData = reshape(modulatedDataSymbol_rotated,[],N_blks).'; % reshape into blocks of 512 samples
            STBC_encoderInputGa64 = repmat(Ga64_rotated,N_blks,1); % Ga64 repeated N_blks-times
            % concatenate STBC encoderInput
%             STBC_encoderInputBlock = [STBC_encoderInputData, STBC_encoderInputGa64]; % append rotated Ga64 to data blocks
            STBC_encoderInputBlock = [STBC_encoderInputData]; % append rotated Ga64 to data blocks
            %         STBC_encoderInputStream = reshape(STBC_encoderInputBlock.',1,[]);
            
%             encodedDataSymbolsSTBC =
%             AlDhahir_STBC_Encoder(STBC_encoderInputBlock, wifi_params); %  puvodni verze - problem s ISI (UW)
            [yblock_2x1_A_tmp, yblock_2x1_B_tmp] = AlDhahir_STBC_Encoder2(N_blks, n_symb, n_gi, STBC_encoderInputBlock); % nova verze, jina ramcova struktura, problemy s ISI castecne vyreseny
            % add UW GI
            if wifi_params.general.dataSymbolRotation == 1
                uw1 = Ga64_rotated;
                uw2 = Gb64_rotated;
            else
                uw1 = Ga64;
                uw2 = Gb64;
            end
            % [zeros, c_block1A, UWa, c_block_2A, zeros]
            % [zeros, c_block1B, UWb, c_block_2B, zeros]
            yblock_2x1_A = [repmat([zeros(1,n_gi); uw1],N_blks/2,1), yblock_2x1_A_tmp(:,:)];
            yblock_2x1_B = [repmat([zeros(1,n_gi); uw2],N_blks/2,1), yblock_2x1_B_tmp(:,:)];
            % reshaping back to row vector
            y_2x1_A = reshape(yblock_2x1_A.', 1, []);
            y_2x1_B = reshape(yblock_2x1_B.', 1, []);
            % add the last GI (= zeros)
            y_2x1_A = [y_2x1_A, zeros(1,n_gi)].'; % column vector
            y_2x1_B = [y_2x1_B, zeros(1,n_gi)].'; % column vector
            % concatenate the output            
            outputBlocks = [y_2x1_A, y_2x1_B]; % nSymbols x nTxAnt (column vector/matrix)
            
        case 'MIMO' % =====================================================
            if nTxAnt ~= 2
                error;
            end
            
            %         modulatedDataSymbol_rotated = (1:length(modulatedDataSymbol_rotated)).';
            % blocks division
            modulatedDataSymbol_blocks = reshape(modulatedDataSymbol_rotated,[],N_blks).';
            % guard insertion (types: GI, CP or zeros)
            GIType = wifi_params.general.guardInterval.type;
            if strcmp(GIType,'GI') == 1
                Ga64_rotated_repeated = repmat(Ga64_rotated, N_blks, 1);
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [Ga64_rotated_repeated, modulatedDataSymbol_blocks];
                theLastGI = Ga64_rotated;
            elseif strcmp(GIType,'CP') == 1
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [modulatedDataSymbol_blocks(:,N_SPB-n_gi+1:N_SPB), modulatedDataSymbol_blocks];
                theLastGI = modulatedDataSymbol_blocks(N_blks,N_SPB-n_gi+1:N_SPB);
            else % zeros only
                % concatenate guard and blocks
                modulatedDataSymbolsAndGuard = [zeros(N_blks, n_gi), modulatedDataSymbol_blocks];
                theLastGI = zeros(1,n_gi);
            end
            
            mimoMapping = 'horizontal'; % 'horizontal' or 'vertical', see IEEE 802.11-16/0388-r0, slide 5 (ppt)
            %             mimoMapping = 'vertical';
            N_antPort = 2;
            N_parallelBlks = (N_blks/N_antPort);
            
            % multiple streams of guard and blocks
            if strcmp(mimoMapping, 'horizontal')
                outputBlocksTmpA = reshape(modulatedDataSymbolsAndGuard(1:2:end,:).',1,[]);
                outputBlocksTmpB = reshape(modulatedDataSymbolsAndGuard(2:2:end,:).',1,[]);
            elseif strcmp(mimoMapping, 'vertical')
                outputBlocksTmpA = reshape(modulatedDataSymbolsAndGuard(:,1:2:end).',1,[]);
                outputBlocksTmpB = reshape(modulatedDataSymbolsAndGuard(:,2:2:end).',1,[]);
            else
                error;
            end
            
            outputBlocks = [outputBlocksTmpA; outputBlocksTmpB];
            
            %             plot(outputBlocks(1,:)), hold all, plot(outputBlocks(2,:))
            % add the last guard on the end of the final blocks
            QQ = (1/sqrt(2)); % precoding matrix QQ = (1/sqrt(2))*eye(2)
            %             QQ = 1;
            outputBlocks = QQ*[outputBlocks, repmat(theLastGI,N_antPort,1)].';
            
            
        otherwise
            
            error;
    end
    
    
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
        Header = [];
        HeaderCRC = [];
        
        
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
        
        % AY: some bug is here:
        modulatedHeaderSymbol_rotated_withGI12 = [modulatedHeaderSymbol_rotated_withGI12, modulatedHeaderSymbol_rotated_withGI12(1,1:1024-length(modulatedHeaderSymbol_rotated_withGI12))]; % fix it
        
        % Beamforming training field (BTF) - in future
        BTF = [];
        modulatedBTFSymbol_rotated = BTF; % Header.*exp(1j*pi*modulated_Header_symbol_index_k/2);
        
        if wifi_params.general.sendDataOnly == 0
            % x_s preallocation
            x_s = zeros(nTxAnt, sum(wifi_params.framing.ChipLengths));
            x_s(:,find(wifi_params.framing.ChipMap{1} == 1)) = 1*repmat(modulatedSTFSymbol_rotated,nTxAnt,1);
            x_s(:,find(wifi_params.framing.ChipMap{2} == 1)) = 1*repmat(modulatedCEFSymbol_rotated,nTxAnt,1);
            x_s(:,find(wifi_params.framing.ChipMap{3} == 1)) = 1*repmat(modulatedHeaderSymbol_rotated_withGI12,nTxAnt,1);
            x_s(:,find(wifi_params.framing.ChipMap{4} == 1)) = outputBlocks.';
            x_s(:,find(wifi_params.framing.ChipMap{5} == 1)) = 0*repmat(modulatedBTFSymbol_rotated,nTxAnt,1);
        else % debug, send only data field
            x_s = outputBlocks.';
        end
        
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


end

