function [hEstPadded, HEst] = ad_cir_estimator(y_rx_cef, y_rx_stf, channelObj, wifi_params, MCS)
%ad_cir_estimator Channel impulse response estimator based on Golay COmplementary sequences
% Author: Jiri Milos, DREL FEEC BUT
% Date: 2019/05/25
verze = 1;

n_gi = wifi_params.mapping.N_GI; % samples

% remove -pi/2 rotation
modulated_CEF_symbol_index_k = (0:length(y_rx_cef)-1); % Channel estimation field
modulated_STF_symbol_index_k = (0:length(y_rx_stf)-1); % Short training field

y_rx_cef_derotated = y_rx_cef.*exp(-1j*pi*modulated_CEF_symbol_index_k/2);
y_rx_stf_derotated = y_rx_stf.*exp(-1j*pi*modulated_STF_symbol_index_k/2);

% get STF suffix (last 128 samples, equal to -Ga128 sequence)
y_rx_stf_128_suffix = y_rx_stf_derotated(end-128+1:end);


% simple Golay method to estimate CIR
Ga128 = wifi_params.spreading.Golay_Seq.Ga_128;
Gb128 = wifi_params.spreading.Golay_Seq.Gb_128;

Ga256 = [Ga128, Gb128];
Gb256 = [Ga128, -Gb128];

Ga512 = [Ga256, Gb256];
Gb512 = [Ga256, -Gb256];


% Gu512 = [-Gb128, -Ga128, +Gb128, -Ga128];
% Gv512 = [-Gb128, +Ga128, -Gb128, -Ga128];
% Gv128 = -Gb128;
% CEF = [Gu512, Gv512, Gv128];

y_rx_cef_128_derotated = reshape(y_rx_cef_derotated,[],9).';

switch verze
    case 1 % puvodni (ziskana delka odezvy = 128 vzorku)
        
        R1 = xcorr(y_rx_cef_128_derotated(1,:), -Gb128);
        R2 = xcorr(y_rx_cef_128_derotated(2,:), -Ga128);
        R3 = xcorr(y_rx_cef_128_derotated(3,:), +Gb128);
        R4 = xcorr(y_rx_cef_128_derotated(4,:), -Ga128);
        R5 = xcorr(y_rx_cef_128_derotated(5,:), -Gb128);
        R6 = xcorr(y_rx_cef_128_derotated(6,:), +Ga128);
        R7 = xcorr(y_rx_cef_128_derotated(7,:), -Gb128);
        R8 = xcorr(y_rx_cef_128_derotated(8,:), -Ga128);
        R9 = xcorr(y_rx_cef_128_derotated(9,:), -Gb128);
        
        R12 = R1 + R2;
        R23 = R2 + R3;
        R34 = R3 + R4;
        R45 = R4 + R5;
        R56 = R5 + R6;
        R67 = R6 + R7;
        R78 = R7 + R8;
        R89 = R8 + R9;
        
        Rab = (1/(2*length(Ga128)))*(sum([R12; R23; R34; R45; R56; R67; R78; R89],1)/8);
        hEst = Rab(length(Ga128):end);
        
        % figure(5)
        % subplot(241)
        % stem(R12), grid on;
        % subplot(242)
        % stem(R23), grid on;
        % subplot(243)
        % stem(R34), grid on;
        % subplot(244)
        % stem(R45), grid on;
        % subplot(245)
        % stem(R56), grid on;
        % subplot(246)
        % stem(R67), grid on;
        % subplot(247)
        % stem(R78), grid on;
        % subplot(248)
        % stem(R89), grid on;
        %
        % figure(6)
        % stem(Rab)
        % hold on
        % stem(128,0,'rx')
    case 2 % pokus o ziskani prodlouzene verze zretezenim Ga a Gb z CES field, pr. Ga256 = [Ga128, Gb128], Gb256 = [Ga128, -Gb128]
        R1signMap = [-1*ones(1, 128), +1*ones(1, 128)];
        R2signMap = [-1*ones(1, 128), +1*ones(1, 128)];
        R3signMap = [+1*ones(1, 128), -1*ones(1, 128)];
        R4signMap = [-1*ones(1, 128), +1*ones(1, 128)];
        
        R1 = xcorr([y_rx_cef_128_derotated(2,:), y_rx_cef_128_derotated(3,:)].*R1signMap, Ga256.*R1signMap);
        R2 = xcorr([y_rx_cef_128_derotated(4,:), y_rx_cef_128_derotated(5,:)].*R2signMap, Gb256.*R2signMap);
        R3 = xcorr([y_rx_cef_128_derotated(6,:), y_rx_cef_128_derotated(7,:)].*R3signMap, Ga256.*R3signMap);
        R4 = xcorr([y_rx_cef_128_derotated(8,:), y_rx_cef_128_derotated(9,:)].*R4signMap, Gb256.*R4signMap);
        
        R12 = R1 + R2;
        R34 = R3 + R4;
        
        Rab = (1/(2*length(Ga256)))*(sum([R12; R34],1)/2);
        hEst = Rab(length(Ga256):end);
    case 3
        
        %         R1 = xcorr([y_rx_stf_128_suffix,         y_rx_cef_128_derotated(1,:)], [+Ga128, -Gb128]); % corr seq. A
        %         R2 = xcorr([y_rx_cef_128_derotated(1,:), y_rx_cef_128_derotated(2,:)], [-Gb128, -Ga128]); % B
        %         R3 = xcorr([y_rx_cef_128_derotated(2,:), y_rx_cef_128_derotated(3,:)], [-Ga128, +Gb128]); % A
        %         R4 = xcorr([y_rx_cef_128_derotated(3,:), y_rx_cef_128_derotated(4,:)], [+Gb128, -Ga128]); % B
        %         R5 = xcorr([y_rx_cef_128_derotated(4,:), y_rx_cef_128_derotated(5,:)], [-Ga128, -Gb128]); % A
        %         R6 = xcorr([y_rx_cef_128_derotated(5,:), y_rx_cef_128_derotated(6,:)], [-Gb128, +Ga128]); % B
        %         R7 = xcorr([y_rx_cef_128_derotated(6,:), y_rx_cef_128_derotated(7,:)], [+Ga128, -Gb128]); % A
        %         R8 = xcorr([y_rx_cef_128_derotated(7,:), y_rx_cef_128_derotated(8,:)], [-Gb128, -Ga128]); % B
        %         R9 = xcorr([y_rx_cef_128_derotated(8,:), y_rx_cef_128_derotated(9,:)], [-Ga128, -Gb128]); % A
        %
        R1 = xcorr([y_rx_stf_128_suffix,         y_rx_cef_128_derotated(1,:)], [+Ga128, -Gb128]); % corr seq. A
        R2 = xcorr([y_rx_cef_128_derotated(2,:), y_rx_cef_128_derotated(3,:)], [-Ga128, +Gb128]); % B
        R3 = xcorr([y_rx_cef_128_derotated(4,:), y_rx_cef_128_derotated(5,:)], [-Ga128, -Gb128]); % A
        R4 = xcorr([y_rx_cef_128_derotated(6,:), y_rx_cef_128_derotated(7,:)], [+Ga128, -Gb128]); % B
        %         R8 = xcorr([y_rx_cef_128_derotated(7,:), y_rx_cef_128_derotated(8,:)], [-Gb128, -Ga128]); % B
        
        
        % get correlations
        %         R12 = R1 + R2; % corr seq. A + B
        %         R23 = R2 + R3; % B + A
        %         R34 = R3 + R4; % A + B
        %         R45 = R4 + R5; % B + A
        %         R56 = R5 + R6; % A + B
        %         R67 = R6 + R7; % B + A
        %         R78 = R7 + R8; % A + B
        %         R89 = R8 + R9; % B + A
        %
        R12 = R1 + R2; % corr seq. A + B
        R34 = R3 + R4; % A + B
        
        %         R89 = R8 + R9; % B + A
        
        
        
        % Rab = (1/(2*length(Ga256)))*(sum([R12; R23; R34; R45; R56; R67; R78; R89],1)/8);
        Rab = (1/(2*length(Ga256)))*(sum([R12; R34],1)/2);
        hEst = Rab(length(Ga256):end);
        
    otherwise
        error('Dosud nedefinovany typ odhadu');
end


% get channel CIR and CTF
switch MCS.phy_type
    case 'Ctrl' % Control -----------------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            hEstPadded = [channelObj.h, zeros(1, 32-length(channelObj.h))];
            HEst = channelObj.H(1:32);
        else
            hEstPadded = hEst(1:32);
            HEst = fft(hEstPadded);
        end
        
    case 'SC' % SC ------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            hEstPadded = [channelObj.h, zeros(1, 448-length(channelObj.h))];
            HEst = channelObj.H(1:end-n_gi);
        else
            hEstPadded = [[hEst(1:120),zeros(1,length(hEst)-120)], zeros(1,n_fft-n_gi-length(hEst))];
            HEst = fft(hEstPadded);
        end
        
    case 'OFDM' % OFDM ------------------------------------------------------
        if strcmp(wifi_params.general.channel,'awgn') == 1
            hEstPadded = [channelObj.h, zeros(1, n_fft-length(channelObj.h))];
            HEst = channelObj.H(wifi_params.mapping.ir_data,:); % 336 symbol positions
        else
            hEstPadded = [[hEst(1:120),zeros(1,length(hEst)-120)], zeros(1,n_fft-length(hEst))];
            %             H_est = fft(est_h_padded);
            %             H_est = H_est(wifi_params.mapping.ir_data);
            %             H_est = H_est(1:336);
            HEst = channelObj.H(wifi_params.mapping.ir_data);
        end
        
    case 'LPSC'
        %         H_est % ma byt bez 8*7 symbolu (mozna bude potreba i jina uprava)
        if strcmp(wifi_params.general.channel,'awgn') == 1
            hEstPadded = [channelObj.h, zeros(1, 56*7-length(channelObj.h))];
            HEst = channelObj.H(1:56*7);
        else
            hEstPadded = [[hEst(1:120),zeros(1,length(hEst)-120)], zeros(1,56*7-128)];
            HEst = fft(hEstPadded);
        end
        
    otherwise
        error('Unsupported physical layer type');
end


end

