function [modulatedDataSymbol_rotated] = ayModulator(modulatorIn, wifi_params, MCS)

switch MCS.phy_type
    case 'Ctrl' % modulation + spreading
        % nondifferential stream, Sec. 20.4.3.3.4
        s_k = 2*modulatorIn-1;
        
        d_k = zeros(size(s_k));
        for i_k = 1:length(s_k)
            if i_k == 1
                d_km1 = 1;
            else
                d_km1 = d_k(i_k-1);
            end
            d_k(i_k,1) = s_k(i_k)*d_km1;
        end
        
        % spreading, Sec. 20.4.3.3.5
        n = (0:(length(d_k)*32)-1).';
        n_mod32 = mod(n, 32);
        N_Spr_Blks = length(n)/32;
        Ga32 = wifi_params.spreading.Golay_Seq.Ga_32;
        Ga32_spr = Ga32(n_mod32+1).';
        i_n_flr_div32 = floor(n/32);
        d_k_spr = d_k(i_n_flr_div32+1);
        modulatedDataSymbol = Ga32_spr.*d_k_spr;
        modulatedDataSymbol_rotated = modulatedDataSymbol.*exp(1j*pi*n/2);
        
    case 'SC' % modulation
        hMod = wifi_params.modulation(log2(MCS.M)).hMod;
                
        if wifi_params.general.useNUC == 0  % common modulation type
            
            modulatedDataSymbol = step(hMod,modulatorIn);
            
        else % Non-Uniform Constellation, TBD
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
            hNUCMod = wifi_params.NUC(idxNUC).hNUCMod;
            
            SymbolAlphabet = wifi_params.NUC(idxNUC).SymbolAlphabet.';
%             SymbAlphabet = wifi_params.NUC(idxNUC).SymbolAlphabet;
%             %             TBD: Interleaver!!!
             inputNucModulatorBin = reshape(modulatorIn,6,[]).';
             inputNucModulatorDec = bi2de(inputNucModulatorBin,'left-msb');
%             
            modNorm = 1;
%             % common modulation
%             modulatedDataSymbol = step(hMod,modulatorIn.')*modNorm;
%             modulatedDataSymbol = (SymbolAlphabet(inputNucModulatorDec+1))*modNorm;
%             inputNucModulatorDec = reshape(de2bi(modulatorIn.',log2(length(SymbolAlphabet)),'left-msb').',1,[]).';
            modulatedDataSymbol = step(hNUCMod, inputNucModulatorDec);
%             
            

%     
%             hQAMMod = comm.GeneralQAMModulator(...
%                 'Constellation',SymbolAlphabet);
%             modulatedDataSymbol = step(hQAMMod, inputNucModulatorDec);

        end
        
        % add pi/2 rotation
        if wifi_params.general.dataSymbolRotation == 0
            modulatedDataSymbol_rotated = modulatedDataSymbol;
        else
            modulated_STF_symbol_index_k = (0:length(modulatedDataSymbol)-1).';
            modulatedDataSymbol_rotated = modulatedDataSymbol.*exp(1j*pi*modulated_STF_symbol_index_k/2);
        end
        
        
        
    case 'OFDM'
        N_mod_groups = numel(modulatorIn)/MCS.N_cbps;
        % broke the input stream into groups of N_cbps bits
        modulatorInputGroups = reshape(modulatorIn,[],N_mod_groups).';
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

end

