function chanObj = channel(txObj, wifi_params, i_mcs, nTxAnt, nRxAnt, SNRdB, i_snr)
% Channel function for creating awgn or fading channel matrix
%
% Authors:	Jiri Milos, DREL FEEC BUT, 2018--2020
%			Jiri Blumenstein, DREL FEEC BUT, 2019--2020
%           Radim Zedka, DREL FEEC BUT, 2020

%% Auxiliary parameters
MCS = wifi_params.MCS(i_mcs);

%% CHANNEL MATRIX
n_fft = wifi_params.mapping.n_fft;
n_gi = wifi_params.mapping.N_GI;
n_symb = n_fft-n_gi;

Fs = wifi_params.Fs;
% % H = zeros(n_fft,1,nTxAnt,nRxAnt);
% H = zeros(n_fft, nTxAnt, nRxAnt);

% if strcmp(wifi_params.general.PHYlayer,'OFDM')
%     Nsym = wifi_params.mapping.n_data_symbols_per_frame;
% end

if strcmp(wifi_params.channel.type,'awgn')
    % generate channel transfer function (CTF)
    H = ones(n_fft, nTxAnt, nRxAnt);
    % channel impulse response
    h = ifft(H);
    
    % save parameters to chanMod
    chParams.normh = 1;
    chParams.h_length = 1;
    chParams.nTX = nTxAnt;
    chParams.nRX = nRxAnt;
    chParams.Lch = [];
    chParams.Ntap = [];
    % aux parameters
    chGenie.PDB_dB = [];
    chGenie.PDB_lin = [];
    chGenie.delay = [];
    chGenie.powlin_for_h_est = abs(h);
    chGenie.h = h;
    chGenie.H = H;
    
elseif strcmp(wifi_params.channel.type,'fad') % fading
    
    if nTxAnt > 1
        error
    end
    
    %     Rayleigh (flat)
    Ntap = 64;
    %     Ntap = 70;
    htmp = (randn(Ntap,1) + 1j*randn(Ntap,1))/sqrt(2);
    h = [htmp; zeros(n_fft-length(htmp),1)];
    normalize_factor_h = max(abs(h));
    h = h/normalize_factor_h;
    
    
    % for perfect channel knowledge
    H = fft(h);
    
    
    
    
    % % according to taps
    %     error('Undefined yet');
    %     powerdB = wifi_params.channel.powerdB;
    %     powerlin = 10.^(powerdB/10);
    %
    %     normalize_factor_h = sqrt(sum(10.^powerdB/10));
    %
    %     delayTime = wifi_params.channel.delayTime; % channel delay sample
    %     delay = unique(round(delayTime*Fs));
    %
    %     Ntap = length(powerdB); % channel tap number
    %     Lch = delay(end)+1;
    %
    %     g_rand = (randn(1,Ntap)+1j*randn(1,Ntap))/sqrt(2);
    %
    %     h = zeros(1,Lch);
    %
    %     for i_tap = 1:length(delay)
    %         h(1, delay(i_tap)+1) = sqrt(powerlin(i_tap))*g_rand(1, i_tap);
    %     end
    %
    %     h = h/normalize_factor_h;
    %     h_length = length(h);
    
    %     % for perfect channel knowledge
    %     H = fft([h, zeros(1, n_fft-h_length-n_gi)]);
    %     H2 = fft([h, zeros(1, n_fft-h_length)]);
    %
    %     chanMod.H = H; % only block fading
    %     chanMod.h = h;
    %     chanMod.normh = normalize_factor_h;
    %     chanMod.h_length = h_length;
    %     chanMod.nTX = nTxAnt;
    %     chanMod.nRX = nRxAnt;
    %     chanMod.Lch = Lch;
    %     chanMod.Ntap = Ntap;
    %     chanMod.genie.PDB_dB = powerdB;
    %     chanMod.genie.PDB_lin = powerlin;
    %     chanMod.genie.delay = delay;
    %     chanMod.genie.powlin_for_h_est = abs(h);
    %     chanMod.genie.H2 = H;
    %
    
    % save parameters to chanMod
    chParams.normh = [];
    chParams.h_length = length(h);
    chParams.nTX = nTxAnt;
    chParams.nRX = nRxAnt;
    chParams.Lch = [];
    chParams.Ntap = [];
    % aux parameters
    chGenie.PDB_dB = [];
    chGenie.PDB_lin = [];
    chGenie.delay = [];
    chGenie.powlin_for_h_est = abs(h);
    chGenie.h = h;
    chGenie.H = H;
elseif strcmp(wifi_params.channel.type,'los') % LOS fading
    
    if nTxAnt > 2
        error
    end
    
    maxTxRxAnt = max(nRxAnt, nTxAnt);
    fiMatrix = exp(1j*2*pi*rand(maxTxRxAnt));
    alphadB = -30;
    alphaLin = 10^(alphadB/10);
    alphaMatrix = sqrt(alphaLin)*(ones(maxTxRxAnt)-eye(maxTxRxAnt))+eye(maxTxRxAnt);
    HbaseMatrix = alphaMatrix.*fiMatrix;
    HbaseMatrix = HbaseMatrix(1:nRxAnt, 1:nTxAnt);
    
    Ntap = 1;
    
    % create channel matrix
    for iTx = 1:nTxAnt
        for iRx = 1:nRxAnt
            
            h_tmp = HbaseMatrix(iRx, iTx)/sqrt(2);
            h_tmp_max = max(abs(h_tmp));
            h_tmp_norm = h_tmp./h_tmp_max; % CIR magnitude normalization
            
            % Align to the strongest tap:
            %             [~, idx_h_max] = max(abs(h_tmp_norm));
            %             h_tmp2 = h_tmp_norm(idx_h_max:end);
            %             h2 = (h_tmp2*ones(n_fft,1))/sqrt(sum(abs(h_tmp2).^2));
            h2 = (h_tmp_norm*ones(n_fft,1));
            % Normalization coefficient given by sqrt(sum(abs(h_11)^2))
            % Now I am sure that the channel doesn't increase the signal power.
            
            h(:,iTx,iRx) = h2;
            H(:,iTx,iRx) = fft(h2);
            
        end
    end
    
    
    % save parameters to chanMod
    chParams.normh = [];
    chParams.h_length = length(h);
    chParams.nTX = nTxAnt;
    chParams.nRX = nRxAnt;
    chParams.Lch = [];
    chParams.Ntap = [];
    % aux parameters
    chGenie.PDB_dB = [];
    chGenie.PDB_lin = [];
    chGenie.delay = [];
    chGenie.powlin_for_h_est = abs(h);
    chGenie.h = h;
    chGenie.H = H;
    
    
    
elseif strcmp(wifi_params.channel.type,'fad_meas') % fading - measured channel
    
    % load measured CTF
    N = 1001; % Number of CIR (CTF, S21) samples
    f_span = 1e10; % Frequency span of the S21 measurements
    df = f_span/(N-1); % freq resolution [Hz]
    f1 = 0;
    %     f2 = 2e9; % f_max
    f2 = 2.176e9; % f_max
    n_idx = 1:f2/df;
    f_select = f1:df:f2;
    load('.\measured_channels\CTF.mat'); % creates CTF_mat variable
    
    % Assuming the frequency range is <-F/2, F/2> and not <0, F>
    % % % v 1 spise ne
    %     h = zeros(n_fft-n_gi, nTxAnt, nRxAnt);
    %     H = ones(n_fft-n_gi, nTxAnt, nRxAnt);
    
    % % v2 nadejnejsi test prijimca
    %     h = zeros(n_fft, nTxAnt, nRxAnt);
    %     H = ones(n_fft, nTxAnt, nRxAnt);
    %
    % SISO Channel Transfer Function (S21 parameter):
    % create channel matrix
    for iTx = 1:nTxAnt
        for iRx = 1:nRxAnt
            S21tmp = squeeze(CTF_mat(iTx, iRx, n_idx));
            
            h_tmp = ifft(fftshift(S21tmp));
            h_tmp_max = max(abs(h_tmp));
            h_tmp_norm = h_tmp./h_tmp_max; % CIR magnitude normalization
            
            % Align to the strongest tap:
            [~, idx_h_max] = max(abs(h_tmp_norm));
            h_tmp2 = h_tmp_norm(idx_h_max:end);
            %             h2 = [h_tmp2(1:n_gi,1); zeros(n_symb, 1)]/sqrt(sum(abs(h_tmp2).^2));
            
            h2 = [h_tmp2(:,1); zeros(n_fft-length(h_tmp2), 1)]/sqrt(sum(abs(h_tmp2).^2));
            % Normalization coefficient given by sqrt(sum(abs(h_11)^2))
            % Now I am sure that the channel doesn't increase the signal power.
            
            h(:,iTx,iRx) = h2;
            H(:,iTx,iRx) = fft(h2);
            
            
        end
    end
    
    
    
    % save parameters to chanMod
    chParams.normh = [];
    chParams.h_length = [];
    chParams.nTX = nTxAnt;
    chParams.nRX = nRxAnt;
    chParams.Lch = [];
    chParams.Ntap = [];
    % aux parameters
    chGenie.PDB_dB = [];
    chGenie.PDB_lin = [];
    chGenie.delay = [];
    chGenie.powlin_for_h_est = abs(h2);
    chGenie.h = h;
    chGenie.H = H;
    
else
    error('Unsupported channel type (see channel.m)');
    
end

chParams.type = wifi_params.channel.type;
if (nTxAnt == 1) && (nRxAnt == 1) % simple Ant.  config. division
    chParams.confguration = 'SISO';
elseif (nTxAnt == 1) && (nRxAnt > 1)
    chParams.confguration = 'SIMO';
elseif (nTxAnt > 1) && (nRxAnt == 1)
    chParams.confguration = 'MISO';
else
    chParams.confguration = 'MIMO';
end
%% CHANNEL OUTPUT


% transmitter output
x_s = txObj.output.x_s.'; % column vector


var_x_s = var(x_s); % original signal variance
nSymbols = length(x_s(:,1));

% SIMO and MIMO TxD cases - duplicate tx signal for each receive antenna
if (nTxAnt >= 1) && (nRxAnt > 1)
    x_s = repmat(x_s, 1, 1, nRxAnt);
end


% prealloc
x_sh_rxAnt = zeros(nSymbols,1,nRxAnt);

% set correlation matrix (MIMO only, TXcorrMatrix = RXcorrMatrix)
if strcmp(wifi_params.antConfig.mode,'MIMO') == 1
    antCorrFactor = 0.1;
    hCorrMatrix = eye(max([nRxAnt, nTxAnt]))+((ones(max([nRxAnt, nTxAnt]))-eye(max([nRxAnt, nTxAnt])))*antCorrFactor);
else
    hCorrMatrix = ones(nRxAnt, nTxAnt); % TxD, etc.
end


%% Convolution part
for kk = 1:nRxAnt
    % select h and x_s for correct RxAntenna
    h_RxAnt = h(:,:,kk); % impulse response(s) for kk-th antenna
    x_s_RxAnt = x_s(:,:,kk);
    
    % prealloc variable per transmitting anntenna
    x_sh_txAnt = zeros(nSymbols,nTxAnt);
    
    for jj = 1:nTxAnt
        % select h and x_s for correct TxAntenna
        if nTxAnt > 1
            h_txAnt = h_RxAnt(:,jj)/sqrt(2);
        else
            h_txAnt = h_RxAnt(:,jj);
        end
        x_s_txAnt = x_s_RxAnt(:,jj);
        % convolution
        x_sh_filt = conv(x_s_txAnt, hCorrMatrix(kk,jj)*h_txAnt); % impact of correlation matrix added
        x_sh = x_sh_filt(1:nSymbols); % cut off the end of convolution result
        x_sh_txAnt(:,jj) = x_sh;
    end
    % sum convolved signal from several transmitted anttenas for each
    % received antenna
    x_sh_rxAnt(:,1,kk) = sum(x_sh_txAnt,2);
    
end


% channel output

% x_sh_filt = zeros(1, length(txObj.output.x_s)+channelObj.h_length-1);
% x_sh_filt(1,1:end) = conv(txObj.output.x_s,channelObj.h); % puvodni
% channelObj.x_h = x_sh_filt(1:length(txObj.output.x_s)); % useknuti konce

% add AWGN noise to each receiving antenna
% channelObj.SNR = SNR(i_snr);
chGenie.SNR = SNRdB; % save SNR in dB to chanMod structure

% n = 10^(-SNRdB(i_snr)/20)*((randn(size(channelObj.x_h))+1j*(randn(size(channelObj.x_h))))/sqrt(2));
n = 10^(-SNRdB/20)*((randn(size(x_sh_rxAnt))+1j*(randn(size(x_sh_rxAnt))))/sqrt(2));

% if strcmp(MCS.phy_type,'OFDM') == 1
%     v = (1/sqrt(wifi_params.mapping.n_tot.nonzero)) * n;
% else
%     v = sqrt(wifi_params.mapping.n_fft/wifi_params.mapping.n_tot) * n;
% end
v = n;
% noise variance
chGenie.var_v = var(v);
chGenie.var_x_s = var_x_s;

%% PA impairments
usePAimp = 0;
if usePAimp == 1
    if nRxAnt > 1 || nTxAnt > 1
        error('Not implemented for MIMO yet')
    end
    avgPow = 1e0;
    minD = avgPow2MinD(avgPow, 64);
    hMemlesNonlin = comm.MemorylessNonlinearity('Method','Rapp model'); % Rapp model
    x_sh_rxAnt = sqrt(1/(minD/2)).*hMemlesNonlin(x_sh_rxAnt);
end

usePhaNoise = 0;
if usePhaNoise == 1
    if nRxAnt > 1 || nTxAnt > 1
        error('Not implemented for MIMO yet')
    end
    hPhzNoise = comm.PhaseNoise('Level',[-90 -130],'FrequencyOffset',[1e6 100e6],'SampleRate',wifi_params.Fs);
    x_sh_rxAnt = hPhzNoise(x_sh_rxAnt);
end

%% Save to output structure
chanObj.output.x_hn = x_sh_rxAnt+v;
chanObj.genie = chGenie;
chanObj.params = chParams;

end

function minD = avgPow2MinD(avgPow, M)
    % Average power to minimum distance    
    nBits = log2(M);
    if (mod(nBits,2)==0)
        % Square QAM
        sf = (M - 1)/6;
    else
        % Cross QAM
        if (nBits > 4)
            sf = ((31 * M / 32) - 1) / 6;
        else
            sf = ((5 * M / 4) - 1) / 6;
        end
    end
    minD = sqrt(avgPow/sf);
end
