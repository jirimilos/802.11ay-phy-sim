function dataDemodulatorLLR = ayDemodulator(rx_modulatedDataSymbol_blocks, chanObj, HEst, wifi_params, MCS, noiseVar)

% used demodulator block (Approximate log-likelihood ratio)

cbpb_index = 2; % normal length, TEST
n_fft = wifi_params.mapping.n_fft;
n_gi = wifi_params.mapping.N_GI;
n_symb = n_fft-n_gi;


nTxAnt = wifi_params.antConfig.nTxAnt;
nRxAnt = wifi_params.antConfig.nRxAnt;

% processing per blocks in SC, OFDM and Ctrl, and in signle row for LPSC
switch MCS.phy_type
    case {'Ctrl'} % Control -----------------------------------------------------------------
        hDemod = wifi_params.modulation(log2(MCS.M)).hDemodCtrl;
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb);
        
        Q = 32; % number of chips
        n_fft_spr = 2^(nextpow2(Q)+1);
        p = (n_fft_spr-Q)/2;
        
        % auxiliary Golay sequences
        Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
        %     n = (0:y_rx_data_length-1).';
        %     rx_modulatedDataSymbol_derotated = y_rx_data.*repmat(exp(-1j*pi*n/2),1,1,nRxAnt);
        
        
        n_windows = size(rx_modulatedDataSymbol_blocks,1);
        
        hEst = ifft(HEst);
        HEst_spr = fft(hEst(1:64));
        
        % overlap-cut method
        for iw = 1:n_windows
            slidingWindowInput = fft(rx_modulatedDataSymbol_blocks(iw,:));
            % input equalized symbols
            rxEqualizedTmp = equalizerSISO(slidingWindowInput, HEst_spr.', noiseVar, 'MMSE');
            rxEqualizedSumRotated = ifft(rxEqualizedTmp(1,p+(1:32))); %
            rxEqualizedSum = rxEqualizedSumRotated; %.*exp(-1j*pi*(0:Q-1)/2);
            
            rxEqualized(1,iw) = sum(rxEqualizedSum.*Ga32)/Q;
        end
        
%         konstelace(rxEqualized)
        % auxiliary Golay sequences
        Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
        % Despreading + equalization
        % prealloc
        N_spreadedSymbols = size(rx_modulatedDataSymbol_blocks, 1);
        rxDespreaded = zeros(N_spreadedSymbols,1);
        for iSprSymb = 1:N_spreadedSymbols
            % equalizer input word
            FDE_dataInput = fft(rx_modulatedDataSymbol_blocks(iSprSymb,:)).';
            % input equalized symbols
            rxEqualized = equalizerSISO(FDE_dataInput, HEst(1:32), noiseVar, 'MMSE');
            rxDespreadedTmp = ifft(rxEqualized).*Ga32.';
            rxDespreaded(iSprSymb,1) = sum(rxDespreadedTmp,1)/32;
            
        end
        
        % differential demapping
        s_k = zeros(size(rxDespreaded));
        for i_k = 1:length(s_k)
            if i_k == 1
                d_km1 = 1;
            else
                d_km1 = (rxDespreaded(i_k-1));
            end
            s_k(i_k,1) = rxDespreaded(i_k)*d_km1;
        end
        
        dataDemodulatorLLR = step(hDemod, s_k, noiseVar);
        
        
        %         dataDemodulatorLLR_HD = (dataDemodulatorLLR > 0);
        
    case {'SC'}
        
        hDemod = wifi_params.modulation(log2(MCS.M)).hDemod;
        hSphDecoder = wifi_params.modulation(log2(MCS.M)).hDemodSSD; % soft-sphere decoder
        
        % number of codewords
        N_cw = wifi_params.coding.N_cw;
        % single codeword length
        L_ldpc = wifi_params.coding.L_ldpc;
        % number of symbol blocks
        N_blks = ceil((N_cw*L_ldpc)/MCS.N_cbpb(cbpb_index));
        
        dataDemodulatorLLR_blk = zeros(log2(MCS.M)*448, N_blks);
        
        if nTxAnt == 1
            
            for iBlk = 1:N_blks
                %                 idxBlk = ((iBlk-1)*448)+(1:448);
                % equalizer input word
                FDE_dataInput = fft(rx_modulatedDataSymbol_blocks(iBlk,:,:));
                
                if wifi_params.general.useNUC == 1 % Non-Uniform Constellation
                    switch MCS.CR
                        case 1/2
                            idxNUC = 1;
                        case 5/8
                            idxNUC = 2;
                        case 3/4
                            idxNUC = 3;
                        case 13/16
                            idxNUC = 4;
                        case 7/8 % test: experimental modulation mapping
                            idxNUC = 5;
                    end
                    hNUCDemod = wifi_params.NUC(idxNUC).hNUCDemod;
                    %                     modNorm = 1/sqrt(2);
                    modNorm = 1;
                    SymbolAlphabet = modNorm*wifi_params.NUC(idxNUC).SymbolAlphabet.';
                    bittable = (wifi_params.NUC(idxNUC).bittable);
%                         bittable = double(fliplr(wifi_params.NUC(idxNUC).bittable)).';

                                

                end
                
                
                % input equalized symbols
                switch chanObj.params.confguration
                    case 'SISO'
                        if strcmp(chanObj.params.type,'awgn')
                            rxEqualizedTDRot = (equalizerSISO(FDE_dataInput, HEst(:,:,:).', noiseVar, 'ZF'));
                        else                            
                            rxEqualizedTDRot = (equalizerSISO(FDE_dataInput, HEst(:,:,:).', noiseVar, 'MMSE'));
                        end
                        %                         rxEqualized = ifft(equalizerSISO(FDE_dataInput, HEst(:,:,:), noiseVar, 'ZF'));
                    case 'SIMO'
                        rxEqualizedTDRot = (equalizerSIMO(squeeze(FDE_dataInput), squeeze(HEst(:,:,:)), noiseVar, 'MMSE')).';
                end
                
                rxEqualizedRot = ifft(rxEqualizedTDRot);
                
                % data part derotation
                if wifi_params.general.dataSymbolRotation == 1
                    modulated_symbol_index_k = (0:n_fft-n_gi-1);
                    rxEqualizedDerot = rxEqualizedRot(1:n_symb).*exp(-1j*pi*modulated_symbol_index_k/2);
                else
                    rxEqualizedDerot = rxEqualizedRot(1:n_symb);
                end
                
                % GI part derotation
                if wifi_params.general.dataSymbolRotation == 1
                    k_Ga64 = 0:n_gi-1;
                    GIDerot = rxEqualizedRot(n_symb+1:end).*exp(-1j*pi*k_Ga64/2);
                else
                    GIDerot = rxEqualizedRot(n_symb+1:end);
                end
                
                %                 if iBlk == N_blks
                %                     konstelace(rxEqualizedDerot, GIDerot), drawnow
                %                     pause(2)
                %                     close all
                %                 end
                %
                if wifi_params.general.useNUC == 1 % Non-Uniform Constellation
%                     dataDemodulatorLLR_blk(:, iBlk) = -step(hNUCDemod, rxEqualizedDerot.', noiseVar);
                    dataDemodulatorLLR_blk(:, iBlk) = step(hNUCDemod, rxEqualizedDerot.', noiseVar);
                else
                    dataDemodulatorLLR_blk(:, iBlk) = step(hDemod, rxEqualizedDerot.', noiseVar);
                end
                %             %
                %             CSM = hDemod.CustomSymbolMapping;
                %             if (MCS.M == 2)
                %                 bittable = (logical(de2bi(0:length(symbolAlphabet)-1).'));
                %             else
                %                 bittable = flipud(logical(de2bi(CSM).'));
                %             end
                %
                %             FDE_dataInput = (rx_modulatedDataSymbol_blocks_derotated(i_eq,:));
                %
                
                
                %             for jj = 1:length(FDE_dataInput)
                %                 puvodni   = 0; % puvodni demodulator (pouze pro 16QAM - overeni)
                %                 if puvodni == 1
                %                     rx_user_sym_tmp = fft(FDE_dataInput(jj));
                %                     H_est_tmp = H_est(jj);
                %
                %                     rx_layer_x(i_eq, jj) = ifft((H_est_tmp.^(-1)).*rx_user_sym_tmp); % ZF
                %                     %
                %                     [Q, R] = qr(H_est_tmp);
                %
                %                     demodulatorOutSSD_blk(1, ((jj-1)*log2(MCS.M))+1:jj*log2(MCS.M)) = (LTE_softsphere(rx_layer_x(i_eq, jj), (rx_user_sym_tmp), Q, R, SymbolAlphabet, bittable, 1, log2(MCS.M))).';
                %                 else
                %                     symbolAlphabet = hDemod.constellation.'; % definice vsech moznych symbolu dane modulace (vstupni parametr SSDekoderu)
                %                     bittable = de2bi(wifi_params.modulation(log2(MCS.M)).customSymbolMapping,'left-msb').'; % definice prislusnych bitovych hodnot parametru symbolAlphabet (vstupni parametr SSDekoderu)
                %                     H_est_tmp = H_est(jj);
                %                     rx_user_sym_tmp = fft(FDE_dataInput(jj));
                %                     rx_layer_x(i_eq, jj) = ifft((H_est_tmp.^(-1)).*rx_user_sym_tmp); % ZF
                %                     modFactor = 1;
                %                     demodulatorOutSSD_blk(1, ((jj-1)*log2(MCS.M))+1:jj*log2(MCS.M)) = -step(hDemod, rx_layer_x(i_eq, jj)*modFactor, noiseVar);
                %
                %                 end
                %
                %
                %             end
                %             demodulatorOutSSD(i_eq,:) = demodulatorOutSSD_blk;
                %%% see: http://www.ursi.org/proceedings/procGA11/ursi/C02-3.pdf % jednoducha implementace MMSE
                
            end
            
            
            dataDemodulatorLLR = reshape(dataDemodulatorLLR_blk,[],1);
            
        else % nTxAnt > 1
            
            if nRxAnt == 1 % 2x1
                % 2x1 system Equalization (STBC decoding):
                Y_2x1_decoded = zeros(N_blks, n_symb);
                
                yblock_rx_2x1 = rx_modulatedDataSymbol_blocks;
                
                hEst = ifft(HEst);
                
                
                for i = 1:2:N_blks-1 % processing for two blocks
%                     if strcmpi(wifi_params.channel.type,'awgn') == 1
%                         %                 [Y_aux] = AlDhahir_STBC_Decoder(448, [yblock_rx_2x1(i,:); yblock_rx_2x1(i+1,:)], [h_00; h_10], 1)*sqrt(2);
%                         [Y_aux] = AlDhahir_STBC_Decoder(...
%                             n_fft,...
%                             [yblock_rx_2x1(i,:); yblock_rx_2x1(i+1,:)],...
%                             [hEst(:,1).'; hEst(:,2).'],...
%                             1, wifi_params)*sqrt(2);
%                         % "Y_aux" is the FFT of two data blocks (2 x SyPD)
%                         %                 Y_2x1_decoded(i,:) = Y_aux(1,:)./(abs(H_00).^2 + abs(H_10).^2 + 1/SNR_lin);
%                         %                 Y_2x1_decoded(i+1,:) = Y_aux(2,:)./(abs(H_00).^2 + abs(H_10).^2 + 1/SNR_lin);
%                         Y_2x1_decodedAtmp = ifft(Y_aux(1,:)./(abs(HEst(1:n_fft,1).').^2 + abs(HEst(1:n_fft,2).').^2 + noiseVar));
%                         Y_2x1_decodedBtmp = ifft(Y_aux(2,:)./(abs(HEst(1:n_fft,1).').^2 + abs(HEst(1:n_fft,2).').^2 + noiseVar));
%                         
%                         Y_2x1_decoded(i,:) = Y_2x1_decodedAtmp(1,1:n_symb);
%                         Y_2x1_decoded(i+1,:) = Y_2x1_decodedBtmp(1,1:n_symb);
%                     else
                        
                        [Y_aux] = AlDhahir_STBC_Decoder2(...
                            n_fft,...
                            [yblock_rx_2x1(i,:); yblock_rx_2x1(i+1,:)],...
                            [hEst(:,1).'; hEst(:,2).'],...
                            1);
                        % linear processing w noise variance (MMSE)
                        Y_2x1_decodedAtmp = (Y_aux(1,:)./(abs(HEst(1:n_fft,1).').^2 + abs(HEst(1:n_fft,2).').^2 + noiseVar));
                        Y_2x1_decodedBtmp = (Y_aux(2,:)./(abs(HEst(1:n_fft,1).').^2 + abs(HEst(1:n_fft,2).').^2 + noiseVar));
                        
                        
                        Y_2x1_decodedAB = ifft([Y_2x1_decodedAtmp; Y_2x1_decodedBtmp].').';
                        %                     konstelace(Y_2x1_decodedAB)
                        
                        Y_2x1_decoded(i,:) = Y_2x1_decodedAB(1,n_gi/2+1:end-(n_gi/2));
                        Y_2x1_decoded(i+1,:) = Y_2x1_decodedAB(2,n_gi/2+1:end-(n_gi/2));
%                     end
                    % allows using received Ga64 to set Decision
                    % Feedback-FDE (not defined yet)
                    Y_2x1_decodedGa64A = Y_2x1_decodedAtmp(1,n_symb+1:end); % shall be redefined
                    Y_2x1_decodedGa64B = Y_2x1_decodedBtmp(1,n_symb+1:end);
                end
                % ifft performed above, within each step of the cycle
                
                % Derotation of Data symbol
                
                if wifi_params.general.dataSymbolRotation == 1
                    modulated_symbol_index_k = (0:numel(Y_2x1_decoded(1,:))-1);
                    modulated_symbol_derotating_factorMatrix = repmat(exp(-1j*pi*modulated_symbol_index_k/2),N_blks,1);
                    % derotate
                    yblock_2x1_decoded_derotated = Y_2x1_decoded.*modulated_symbol_derotating_factorMatrix;
                else
                    yblock_2x1_decoded_derotated = Y_2x1_decoded;
                end
                
                rxEqualizedDerot = reshape(yblock_2x1_decoded_derotated.',1,[]).';
%                 konst = hDemod.constellation;
%                                                     konstelace(rxEqualizedDerot, rxEqualizedDerot*sqrt(2), konst), drawnow
%                                                     pause(2)
%                                                     close all
%                 
                %
                                dataDemodulatorLLR = step(hDemod,rxEqualizedDerot*sqrt(2), noiseVar);
                
%                 dataDemodulatorLLR = -step(hSphDecoder, rxEqualizedDerot*sqrt(2), complex(ones(size(rxEqualizedDerot))));
                
                
            elseif nRxAnt == 2
                if strcmp(wifi_params.antConfig.mode,'TxD') == 1
                    % 2x2 system Equalization (STBC decoding):
                    Y_2x2_decoded = zeros(N_blks, 448);
                    
                    %                 yblock_rx_2x1_GI = reshape(y_rx_2x1(1:Nsmp-NH), SyPD+SyPCP, NBLKpF).';
                    %                 yblock_rx_2x1 = reshape(rx_modulatedDataSymbol_blocks_derotated,448,[]).';
                    %
                    %
                    %                 yblock_rx_2x2_A_GI = reshape(y_rx_2x2_A(1:Nsmp-NH), SyPD+SyPCP, NBLKpF).';
                    %          yblock_rx_2x2_B_GI = reshape(y_rx_2x2_B(1:Nsmp-NH), SyPD+SyPCP, NBLKpF).';
                    
                    yblock_rx_2x2_A = rx_modulatedDataSymbol_blocks(:,:,1);
                    yblock_rx_2x2_B = rx_modulatedDataSymbol_blocks(:,:,2);
                    
                    HEst0 = HEst(:,:,1);
                    HEst1 = HEst(:,:,2);
                    
                    hEst0 = ifft(HEst0).';
                    hEst1 = ifft(HEst1).';
                    
                    %                 hEst0 = hEst0(:,1:448);
                    %                 hEst1 = hEst1(:,1:448);
                    
                    for i = 1:2:N_blks-1
                        %                     [Y_aux] = AlDhahir_STBC_Decoder( SyPD, [yblock_rx_2x2_A(i,:); yblock_rx_2x2_A(i+1,:);yblock_rx_2x2_B(i,:); yblock_rx_2x2_B(i+1,:);], [h_00; h_01; h_10; h_11], 11)*sqrt(2);
                        [Y_aux] = AlDhahir_STBC_Decoder(...
                            n_fft,...
                            [yblock_rx_2x2_A(i,:); yblock_rx_2x2_A(i+1,:);yblock_rx_2x2_B(i,:); yblock_rx_2x2_B(i+1,:);],...
                            [hEst0(1,:); hEst0(2,:); hEst1(1,:); hEst1(2,:)],...
                            11)*sqrt(2);
                        % "Y_aux" is the FFT of two data blocks (2 x SyPD)
                        %                     Y_2x2_decoded(i,:) = Y_aux(1,:)./(abs(HEst0(1:448,1).').^2 + abs(HEst0(1:448,2).').^2 + abs(HEst1(1:448,1).').^2 + abs(HEst1(1:448,2).').^2 + noiseVar(1));
                        %                     Y_2x2_decoded(i+1,:) = Y_aux(2,:)./(abs(HEst0(1:448,1).').^2 + abs(HEst0(1:448,2).').^2 + abs(HEst1(1:448,1).').^2 + abs(HEst1(1:448,2).').^2 + noiseVar(1));
                        %
                        Y_2x2_decodedAtmp = ifft(Y_aux(1,:)./(abs(HEst0(1:n_fft,1).').^2 + abs(HEst0(1:n_fft,2).').^2 + abs(HEst1(1:n_fft,1).').^2 + abs(HEst1(1:n_fft,2).').^2 + noiseVar(1)));
                        Y_2x2_decodedBtmp = ifft(Y_aux(2,:)./(abs(HEst0(1:n_fft,1).').^2 + abs(HEst0(1:n_fft,2).').^2 + abs(HEst1(1:n_fft,1).').^2 + abs(HEst1(1:n_fft,2).').^2 + noiseVar(1)));
                        
                        Y_2x2_decoded(i,:) = Y_2x2_decodedAtmp(1:n_symb);
                        Y_2x2_decoded(i+1,:) = Y_2x2_decodedBtmp(1:n_symb);
                        % allows using received Ga64 to set Decision
                        % Feedback-FDE (not defined yet)
                        Y_2x2_decodedGa64A = Y_2x2_decodedAtmp(1,n_symb+1:end);
                        Y_2x2_decodedGa64B = Y_2x2_decodedBtmp(1,n_symb+1:end);
                    end
                    % ifft performed above, within each step of the cycle
                    
                    % Derotation of Data symbol
                    if wifi_params.general.dataSymbolRotation == 1
                        modulated_symbol_index_k = (0:numel(Y_2x2_decoded(1,:))-1);
                        modulated_symbol_derotating_factorMatrix = repmat(exp(-1j*pi*modulated_symbol_index_k/2),N_blks,1);
                        % derotate
                        yblock_2x2_decoded_derotated = Y_2x2_decoded.*modulated_symbol_derotating_factorMatrix;
                    else
                        yblock_2x2_decoded_derotated = Y_2x2_decoded;
                    end
                    
                    y_2x2_decoded = reshape(yblock_2x2_decoded_derotated.', 1, []).';
                    
                    
                    dataDemodulatorLLR = step(hDemod, y_2x2_decoded, noiseVar);
                    
                elseif strcmp(wifi_params.antConfig.mode,'MIMO') == 1
                    
                    QQ = (1/sqrt(2)); % precoding matrix QQ = (1/sqrt(2))*eye(2)
                    %                     QQ = 1;
                    QQinv = pinv(QQ);
                    % inverse precoding
                    rx_modulatedDataSymbol_blocks = QQinv*rx_modulatedDataSymbol_blocks;
                    
                    % get n-th RX antenna input
                    yblock_rx_2x2_A = rx_modulatedDataSymbol_blocks(:,:,1);
                    yblock_rx_2x2_B = rx_modulatedDataSymbol_blocks(:,:,2);
                    
                    N_antPort = 2;
                    N_parallelBlks = (N_blks/N_antPort);
                    
                    mimoMapping = 'horizontal'; % 'horizontal' or 'vertical', see IEEE 802.11-16/0388-r0, slide 5 (ppt)
                    %             mimoMapping = 'vertical';
                    
                    % decoding
                    yblock_rx_2x2_AB = zeros(N_blks,n_fft);
                    
                    rxEqualizedRotA = zeros(N_parallelBlks,n_fft);
                    rxEqualizedRotB = zeros(N_parallelBlks,n_fft);
                    
                    % equalization
                    for iBlk = 1:N_parallelBlks
                        % part A
                        FDE_dataInputA = fft(yblock_rx_2x2_A(iBlk,:,:)); % input equalized symbols/partA
                        rxEqualizedTDRotA = (equalizerSISO(FDE_dataInputA, HEst(:,1,1).', noiseVar, 'MMSE'));
                        rxEqualizedRotA(iBlk,:) = ifft(rxEqualizedTDRotA);
                        % part B
                        FDE_dataInputB = fft(yblock_rx_2x2_B(iBlk,:,:)); % input equalized symbols/partB
                        rxEqualizedTDRotB = (equalizerSISO(FDE_dataInputB, HEst(:,2,2).', noiseVar, 'MMSE'));
                        rxEqualizedRotB(iBlk,:) = ifft(rxEqualizedTDRotB);
                    end
                    
                    % demultiplex streams of blocks
                    if strcmp(mimoMapping, 'horizontal')
                        yblock_rx_2x2_AB(1:2:end,:) = rxEqualizedRotA;
                        yblock_rx_2x2_AB(2:2:end,:) = rxEqualizedRotB;
                        %                         outputBlocksTmpA = reshape(modulatedDataSymbolsAndGuard(1:2:end,:).',1,[]);
                        %                         outputBlocksTmpB = reshape(modulatedDataSymbolsAndGuard(2:2:end,:).',1,[]);
                    elseif strcmp(mimoMapping, 'vertical')
                        yblock_rx_2x2_AB(:,1:2:end) = rxEqualizedRotA;
                        yblock_rx_2x2_AB(:,2:2:end) = rxEqualizedRotB;
                        %                         outputBlocksTmpA = reshape(modulatedDataSymbolsAndGuard(:,1:2:end).',1,[]);
                        %                         outputBlocksTmpB = reshape(modulatedDataSymbolsAndGuard(:,2:2:end).',1,[]);
                    else
                        error;
                    end
                    
                    % remove guard
                    yblock_rx_2x2_woGuard = yblock_rx_2x2_AB(:,1:n_symb);
                    
                    % derotation
                    if wifi_params.general.dataSymbolRotation == 1
                        modulated_symbol_index_k = (0:numel(yblock_rx_2x2_woGuard(1,:))-1);
                        modulated_symbol_derotating_factorMatrix = repmat(exp(-1j*pi*modulated_symbol_index_k/2),N_blks,1);
                        % derotate
                        yblock_2x2_decoded_derotated = yblock_rx_2x2_woGuard.*modulated_symbol_derotating_factorMatrix;
                    else
                        yblock_2x2_decoded_derotated = yblock_rx_2x2_woGuard;
                    end
                    
                    % demodulation
                    
                    y_2x2_decoded = reshape(sqrt(2)*yblock_2x2_decoded_derotated.', 1, []).';
                    
                    useSphereDecoder = 1;
                    if useSphereDecoder == 0
                        dataDemodulatorLLR = step(hDemod, y_2x2_decoded, noiseVar);
                    else
                        %                         bitTable = de2bi(hDemod.CustomSymbolMapping,log2(hDemod.ModulationOrder),'left-msb');
                        %                         hSphDecoder = comm.SphereDecoder(...
                        %                             'Constellation',constellation(hDemod),...
                        %                             'BitTable',bitTable,'DecisionType','Soft');
                        
                        dataDemodulatorLLR = -step(hSphDecoder, y_2x2_decoded, complex(ones(size(y_2x2_decoded))));
                    end
                    
                    
                end
                
            end
        end
        
    case 'OFDM'
        demodulatorOutSSD = [];
        for i_eq = 1:N_sym
            
            % bittable and symbol alphabet
            symbolAlphabet = [1+1j, -1+1j, -1-1j, 1-1j]/sqrt(2);
            bittable = flipud(logical(de2bi([3 1 0 2]).'));
            % samples of whole OFDM  symbol
            FDE_dataInput = rx_modulatedDataSymbol_blocks(:,i_eq);
            FDE_data_input_tmp = zeros(size(FDE_dataInput));
            % processing within single OFDM symbol
            if (strcmp(MCS.MCS, 'MCS13') == 1) || (strcmp(MCS.MCS, 'MCS14') == 1) % SQPSK is used
                FDE_inp_length = length(FDE_dataInput)/2;
                % prealloc
                demodulatorOutSSD_blk = zeros(1, log2(MCS.M)*FDE_inp_length);
                for jj = 1:FDE_inp_length
                    % first and second SQPSK received symbol
                    rx_user_sym_tmp1 = FDE_dataInput(jj);
                    rx_user_sym_tmp2 = conj(FDE_dataInput(jj+FDE_inp_length));
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
                    demodulatorOutSSD_blk(1,((jj-1)*2)+([1, 2])) = (LTE_softsphere(rx_layer_x, rx_user_sym_tmp, Q, R, symbolAlphabet, bittable, 1, log2(MCS.M)));
                end
                demodulatorOutSSD = [demodulatorOutSSD, demodulatorOutSSD_blk];
            else % QPSK, 16QAM and 64QAM
                FDE_inp_length = length(FDE_dataInput);
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
                        d0R = FDE_dataInput(i_c);
                        d1R = FDE_dataInput(i_c+kR);
                        dR = [d0R; d1R];
                        H_est_tmp = [H_est(i_c); H_est(i_c+kR)];
                        
                        [c_outSSD, ~] = DCM_demod(dR, H_est_tmp, 'QPSK');
                        % save to block
                        demodulatorOutSSD_blk(1, (i_c-1)*4+(1:4)) = c_outSSD;
                    end
                elseif MCS.M == 16
                    for i_c = 1:kR
                        % Static tone pairing demapping
                        d0R = FDE_dataInput(i_c);
                        d1R = FDE_dataInput(i_c+kR);
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
                        dR = FDE_dataInput(i_c);
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
        
        
        
    otherwise
        % error()
        
end

% rx_modulatedDataSymbol_blocks_derotated_equalized = FDE_data_output;
% rx_modulatedDataSymbol_blocks_derotated_equalized = rx_layer_x;






end

