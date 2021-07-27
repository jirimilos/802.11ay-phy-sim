%% Script with main simulation loops (single processing/simulation)
%
% Author:	Jiri Milos, DREL FEEC BUT, 2018
%
% targetBER = 1e-3;
for i_snr = 1:length(SNR) % outer loop - SNR loop
    disp([' SNR = ',num2str(SNR(i_snr)),' dB | ',num2str(i_snr),'/',num2str(length(SNR)),' SNR values']);
    
       
    
    
    for i_loop = 1:N_frames % inner loop - number of frames within the loop
        if mod(i_loop,1000) == 0
            disp(['Part ',num2str(i_loop/1000),'/',num2str(N_frames/1000)]);
        end
        
        %% Transmitter ----------------------------------------------------
        switch wifi_standard
            case '802dot11ad'
                txObj = WIFI_TX_ad(wifi_params,txObj,i_mcs); % perform transmitter operations
            case '802dot11ay'
                txObj = WIFI_TX_ay(wifi_params,txObj,i_mcs); % perform transmitter operations
        end
        %% Channel model, noise and impairments
        % TX Phase Noise
        usePhaNoise = 1;
        if usePhaNoise == 1
            if nRxAnt > 1 || nTxAnt > 1
                error('Not implemented for MIMO yet')
            end
            hPhzNoise = comm.PhaseNoise('Level',[-90 -130],'FrequencyOffset',[1e6 100e6],'SampleRate',wifi_params.Fs);
            tmp1 = txObj.output.x_s;
            tOut = hPhzNoise(tmp1.');
            txObj.output.x_s = tOut.';
        end

        % Channel
        chanObj = channel(txObj, wifi_params, i_mcs, txObj.nTxAnt, rxObj.nRxAnt, SNR(i_snr), i_snr);
        
        % RX Phase Noise
        if usePhaNoise == 1
            if nRxAnt > 1 || nTxAnt > 1
                error('Not implemented for MIMO yet')
            end
            hPhzNoise = comm.PhaseNoise('Level',[-90 -130],'FrequencyOffset',[1e6 100e6],'SampleRate',wifi_params.Fs);
            tmp2 = chanObj.output.x_hn;
            rOut = hPhzNoise(tmp2);
            chanObj.output.x_hn = rOut;
        end
        %% Receiver -------------------------------------------------------
        switch wifi_standard
            case '802dot11ad'
                error('');
%                 rxObj = WIFI_RX_ad(wifi_params,rxObj,channelObj, channelObj.x_hn,i_mcs, SNR(i_snr), txObj, var_v, var_txOutput_x_s);
            case '802dot11ay'
                rxObj = WIFI_RX_ay(wifi_params,rxObj,chanObj, i_mcs, SNR(i_snr), txObj);
        end
        results = result_calculation(results,txObj,rxObj,i_snr,i_loop,N_frames);
        
    end % inner loop
    
end % outer loop