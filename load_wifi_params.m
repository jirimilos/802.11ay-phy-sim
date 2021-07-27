%% Script loads WIFI parameters for simulation and compute dependent parameters
%
% Authors:	Jiri Milos and Ladislav Polak, DREL FEEC BUT, 2018--2019
%
if DEBUG_MODE == true
    warning('on')
    type tasks.txt
else
    warning('off')
end
%% load_wifi_params file for load parameters of simulation
% auxiliary parameters for debugging
useScrambling = true;
sendAllZeros = false;
sendAllOnes = false;
dataSymbolRotation = false;
sendDataOnly = false; % turn off other possible fields
showMapping = 0; % display of subcarriers mapping

% anonymous functions
N_oct2bit = @(x) x*8;
fNcw = @(Length, Lcwd) (1+ceil(N_oct2bit(Length-6)/Lcwd));
fLdpcw = @(Length, Ncw) (ceil(N_oct2bit(Length-6)/(Ncw-1)));
fLdplcw = @(Length, Ncw, Ldpcw) (N_oct2bit(Length-6)-(Ncw-2)*Ldpcw);


% other parameters
wifi_params.general.standard = wifi_standard;
wifi_params.general.channel = channelType;
wifi_params.general.channelKnowledge = channelKnowledge;

wifi_params.general.useScrambling = useScrambling;
wifi_params.general.showMapping = showMapping;

wifi_params.general.sendAllZeros = sendAllZeros;
wifi_params.general.sendAllOnes = sendAllOnes;
wifi_params.general.sendDataOnly = sendDataOnly;

wifi_params.general.dataSymbolRotation = dataSymbolRotation;

wifi_params.general.useNUC = useNUC;
wifi_params.general.use8PSK = use8PSK;

wifi_params.general.guardInterval.length = GuardInterval;
wifi_params.general.guardInterval.type = GuardIntervalType;

wifi_params.general.LDPC_lftMtx = LDPCmatrix;


disp('----------------------------------------------------------------')

if (wifi_params.general.sendAllZeros && wifi_params.general.sendAllOnes) == 1
    disp('Conflict: variables sendAllZeros and sendAllOnes are set to TRUE');
    disp('The transmitter will send all data bits as ones');
    disp('----------------------------------------------------------------')
    wifi_params.general.sendAllZeros = false;
    wifi_params.general.sendAllOnes = true;
end

if (wifi_params.general.sendAllZeros || wifi_params.general.sendAllOnes) == 1
    disp('Data scrambling is disabled');
    disp('----------------------------------------------------------------')
    wifi_params.general.useScrambling = false;
end


%% Colission detection and others:
disp('Simulation settings colission detection...'), pause(0.5);
switch wifi_standard
    case '802dot11ad'
        addpath('.\measured_channels');
        if (strcmp(wifi_params.MCS(i_mcs).MCS, 'MCS0') == 1)
            if LENGTH < 14
                LENGTH = 14;
                disp('> MCS0: The LENGTH parameter must be within the range from 14 to 1023');
                disp('> LENGTH was set to 14 - min. possible value.'), disp(' ');
            elseif LENGTH > 1023
                LENGTH = 1023;
                disp('> MCS0: The LENGTH parameter must be within the range from 14 to 1023');
                disp('> LENGTH was set to 1023 - max. possible value.'), disp(' ');
            end
            
        end
        disp('> AD: No colissions found.')
    case '802dot11ay'
        % select pi/2-8psk field applied -> MCS12, MCS13 use 8PSK
        % select pi/2-64nuc field applied -> MCS17 to MCS21 use 64NUC
        % select DCM pi/2-bpsk field applied -> MCS2 to MCS6 use DCM
    otherwise
        error('> Unsupported type of IEEE 802.11 standard');
end

%% PHY layer - selected MCS throughput
switch wifi_standard
    case '802dot11ad'
        maxPHYThroughput = wifi_params.MCS(i_mcs).DataRate;
        
    case '802dot11ay'
        switch GuardInterval
            case 'short'
                ay_gi_idx = 1;
            case 'normal'
                ay_gi_idx = 2;
            case 'long'
                ay_gi_idx = 3;
        end
        N_CB = wifi_params.MCS(i_mcs).N_CB;
        if strcmp(antMode,'MIMO') == 1
            N_ss = 2;
            
        else
            N_ss = 1;
        end
        wifi_params.MCS(i_mcs).N_ss = N_ss;
        maxPHYThroughput = N_CB*N_ss*wifi_params.MCS(i_mcs).DataRate(ay_gi_idx);
        
        
        
        
    case {'802dot11g', '802dot11ac', '802dot11af', '802dot11ah', '802dot11ax'}
        warning('Check in future (see load_wifi_params.m)');
    otherwise
        error('Undefined IEEE 802.11 standard (see load_wifi_params.m)');
end

%% Create TX and RX object
fcar = 60e9;


TxRxAntTmp = TxRxAnt/10;
nRxAntTmp = rem(TxRxAntTmp,1);

nTxAnt = int16(TxRxAntTmp-nRxAntTmp);
nRxAnt = int16(10*nRxAntTmp);

txObj = transmitter(1,nTxAnt,fcar,wifi_params.MCS); % create transmitter object
rxObj = receiver(1,nRxAnt,fcar,wifi_params.MCS); % create receiver object

clear decType

wifi_params.antConfig.mode = antMode;
wifi_params.antConfig.nTxAnt = nTxAnt;
wifi_params.antConfig.nRxAnt = nRxAnt;

%% channel definition
wifi_params.channel.type = channelType;
fadingType = 'block';
wifi_params.channel.fading = fadingType;
wifi_params.channel.powerdB = [0 -9.7]; % channel path relative power
wifi_params.channel.delayTime = [0 35]*10^(-9); % channel path delays

% clear channelType fadingType




%% PHY layer - timing-related constants
switch wifi_standard
    case '802dot11ad' % =====================================================================================================
        if strcmp( wifi_params.MCS(i_mcs).phy_type, 'OFDM') % OFDM
            % basic parameters
            wifi_params.mapping.bandwidth = 2640e6;
            wifi_params.mapping.n_fft = 512; % FFT size
            wifi_params.mapping.n_data = 336;
            wifi_params.mapping.n_pilot = 16;
            wifi_params.mapping.n_dc = 3;
            wifi_params.mapping.N_GI = 64;
            wifi_params.mapping.i_data = sort(kron([-1, 1],[2:9, 11:29, 31:49, 51:69, 71:89, 91:109, 111:129, 131:149, 151:177]));
            wifi_params.mapping.i_pilots = sort(kron([-1, 1],[10, 30, 50, 70, 90, 110, 130, 150]));
            wifi_params.mapping.i_dc = [-1, 0, 1];
            wifi_params.mapping.i_phase_rotated_subc = [];
            wifi_params.mapping.phase_rotation_factor = 1; % 0°
            wifi_params.mapping.n_segment = 1;
            % dependent parameters
            wifi_params.mapping.n_tot.nonzero = wifi_params.mapping.n_data+wifi_params.mapping.n_pilot;
            wifi_params.mapping.n_tot.all = wifi_params.mapping.n_tot.nonzero + wifi_params.mapping.n_dc;
            wifi_params.mapping.df = wifi_params.mapping.bandwidth/wifi_params.mapping.n_fft;
            wifi_params.mapping.bw_real = wifi_params.mapping.n_tot.all*wifi_params.mapping.df; % real occupied bandwidth
            % OFDM mapping real indices (MATLAB notation)
            wifi_params.mapping.ir_offset = (wifi_params.mapping.n_fft/2)+1;
            wifi_params.mapping.ir_all = 1:wifi_params.mapping.n_fft;
            wifi_params.mapping.ir_data = wifi_params.mapping.i_data+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_pilots = wifi_params.mapping.i_pilots+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_dc = wifi_params.mapping.i_dc+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_nonzero = sort([wifi_params.mapping.ir_data, wifi_params.mapping.ir_pilots, wifi_params.mapping.ir_dc]);
            wifi_params.mapping.ir_gi = [1:wifi_params.mapping.ir_nonzero(1)-1, wifi_params.mapping.ir_nonzero(end)+1:wifi_params.mapping.ir_all(end)];
            % OFDM symbol pilot sequence
            wifi_params.mapping.pilot_seq = [-1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1]; % OFDM symbol pilot sequence, see Sec. 20.5.3.2.5
            wifi_params.mapping.ir_phase_rotated_sc = [];
        else % {Ctrl, SC, LPSC}
            % basic parameters
            wifi_params.mapping.bandwidth = 2640e6;
            wifi_params.mapping.n_fft = 512; % FFT size
            wifi_params.mapping.N_GI = 64;
            wifi_params.mapping.N_SPB = 448;
            wifi_params.mapping.Fchip = 1760e6; % chip rate = 2/3*Fs
            wifi_params.mapping.Tchip = 1/wifi_params.mapping.Fchip; % chip time
            wifi_params.mapping.T_seq = 128*wifi_params.mapping.Tchip;
            wifi_params.mapping.T_STF = 17*wifi_params.mapping.T_seq; % detection sequence duration
            wifi_params.mapping.T_CE = 9*wifi_params.mapping.T_seq; % channel estimation sequence duration
            wifi_params.mapping.T_HEADER = 2*512*wifi_params.mapping.Tchip; % header duration
            wifi_params.mapping.F_CPP = 1760e6; % control mode chip rate
            wifi_params.mapping.T_CPP = 1/wifi_params.mapping.F_CPP; % control mode chip time
            wifi_params.mapping.T_STF_CP = 50*wifi_params.mapping.T_seq; % control mode short training field duration
            wifi_params.mapping.T_CE_CP = 9*wifi_params.mapping.T_seq; % control mode channel estimation field duration
            % dependent parameters
            wifi_params.mapping.n_tot = wifi_params.mapping.N_SPB;
            try
                %             wifi_params.mapping.N_BLKS = []; % see 20.6.3.2.3.3
                wifi_params.mapping.T_data = ((wifi_params.mapping.N_BLKS*512)+64)*wifi_params.mapping.Tchip;
            catch
                warning('N_BLKS undefined yet, see 20.6.3.2.3.3 - calculate after LDPC coding definition')
            end
        end
    case '802dot11ay' % =====================================================================================================
        if strcmp( wifi_params.MCS(i_mcs).phy_type, 'OFDM') % OFDM
            % basic parameters
            wifi_params.mapping.bandwidth = 2640e6;
            wifi_params.mapping.n_fft = 512; % FFT size
            wifi_params.mapping.n_data = 336;
            wifi_params.mapping.n_pilot = 16;
            wifi_params.mapping.n_dc = 3;
            wifi_params.mapping.N_GI = 64;
            wifi_params.mapping.i_data = sort(kron([-1, 1],[2:9, 11:29, 31:49, 51:69, 71:89, 91:109, 111:129, 131:149, 151:177]));
            wifi_params.mapping.i_pilots = sort(kron([-1, 1],[10, 30, 50, 70, 90, 110, 130, 150]));
            wifi_params.mapping.i_dc = [-1, 0, 1];
            wifi_params.mapping.i_phase_rotated_subc = [];
            wifi_params.mapping.phase_rotation_factor = 1; % 0°
            wifi_params.mapping.n_segment = 1;
            % dependent parameters
            wifi_params.mapping.n_tot.nonzero = wifi_params.mapping.n_data+wifi_params.mapping.n_pilot;
            wifi_params.mapping.n_tot.all = wifi_params.mapping.n_tot.nonzero + wifi_params.mapping.n_dc;
            wifi_params.mapping.df = wifi_params.mapping.bandwidth/wifi_params.mapping.n_fft;
            wifi_params.mapping.bw_real = wifi_params.mapping.n_tot.all*wifi_params.mapping.df; % real occupied bandwidth
            % OFDM mapping real indices (MATLAB notation)
            wifi_params.mapping.ir_offset = (wifi_params.mapping.n_fft/2)+1;
            wifi_params.mapping.ir_all = 1:wifi_params.mapping.n_fft;
            wifi_params.mapping.ir_data = wifi_params.mapping.i_data+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_pilots = wifi_params.mapping.i_pilots+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_dc = wifi_params.mapping.i_dc+wifi_params.mapping.ir_offset;
            wifi_params.mapping.ir_nonzero = sort([wifi_params.mapping.ir_data, wifi_params.mapping.ir_pilots, wifi_params.mapping.ir_dc]);
            wifi_params.mapping.ir_gi = [1:wifi_params.mapping.ir_nonzero(1)-1, wifi_params.mapping.ir_nonzero(end)+1:wifi_params.mapping.ir_all(end)];
            % OFDM symbol pilot sequence
            wifi_params.mapping.pilot_seq = [-1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1]; % OFDM symbol pilot sequence, see Sec. 20.5.3.2.5
            wifi_params.mapping.ir_phase_rotated_sc = [];
        else % {Ctrl, SC}
            % basic parameters
            wifi_params.mapping.bandwidth = 2640e6;
            wifi_params.mapping.n_fft = 512; % FFT size
            wifi_params.mapping.N_GI = [32 64 128];
            wifi_params.mapping.N_GI = wifi_params.mapping.N_GI(ay_gi_idx); % select 1 item
            wifi_params.mapping.N_SPB = [480 448 384];
            wifi_params.mapping.N_SPB = wifi_params.mapping.N_SPB(ay_gi_idx); % select 1 item
            wifi_params.mapping.Fchip = 1760e6; % chip rate = 2/3*Fs
            wifi_params.mapping.Tchip = 1/wifi_params.mapping.Fchip; % chip time
            wifi_params.mapping.T_seq = 128*wifi_params.mapping.Tchip;
            wifi_params.mapping.T_STF = 17*wifi_params.mapping.T_seq; % detection sequence duration
            wifi_params.mapping.T_CE = 9*wifi_params.mapping.T_seq; % channel estimation sequence duration
            wifi_params.mapping.T_HEADER = 2*512*wifi_params.mapping.Tchip; % header duration
            wifi_params.mapping.F_CPP = 1760e6; % control mode chip rate
            wifi_params.mapping.T_CPP = 1/wifi_params.mapping.F_CPP; % control mode chip time
            wifi_params.mapping.T_STF_CP = 50*wifi_params.mapping.T_seq; % control mode short training field duration
            wifi_params.mapping.T_CE_CP = 9*wifi_params.mapping.T_seq; % control mode channel estimation field duration
            % dependent parameters
            wifi_params.mapping.n_tot = wifi_params.mapping.N_SPB;
            try
                %             wifi_params.mapping.N_BLKS = []; % see 20.6.3.2.3.3
                wifi_params.mapping.T_data = ((wifi_params.mapping.N_BLKS*512)+64)*wifi_params.mapping.Tchip;
            catch
                warning('N_BLKS undefined yet, see 20.6.3.2.3.3 - calculate after LDPC coding definition')
            end
        end
    otherwise
        error('Unsupported type of IEEE 802.11 standard');
end


wifi_params.general.LENGTH = LENGTH;
switch wifi_standard
    case '802dot11ah'
        wifi_params.mapping.service_length = 8;
    otherwise
        wifi_params.mapping.service_length = 16;
end

%% MODEM common params ====================================================
% BPSK --------------------------------------------------------------------
wifi_params.modulation(1).M = 1; % BPSK
wifi_params.modulation(1).k = 2.^wifi_params.modulation(1).M; % number of constellation points
wifi_params.modulation(1).customSymbolMapping = [1 0];
wifi_params.modulation(1).hMod = comm.PSKModulator(...
    'ModulationOrder',wifi_params.modulation(1).k,...
    'PhaseOffset',pi/4,...
    'BitInput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(1).customSymbolMapping);
wifi_params.modulation(1).hDemod = comm.PSKDemodulator(...
    'ModulationOrder',wifi_params.modulation(1).k,...
    'PhaseOffset',pi/4,...
    'BitOutput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(1).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(1).hDemodSSD = comm.SphereDecoder(...
    'Constellation',constellation(wifi_params.modulation(1).hDemod),...
    'BitTable',de2bi(wifi_params.modulation(1).hDemod.CustomSymbolMapping,log2(wifi_params.modulation(1).hDemod.ModulationOrder),'left-msb'),...
    'DecisionType','Soft');
wifi_params.modulation(1).hDemodCtrl = comm.PSKDemodulator(... % demodulator for Control (Ctrl) PHY
    'ModulationOrder',wifi_params.modulation(1).k,...
    'PhaseOffset',0,...
    'BitOutput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(1).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(1).hMod_factor = 1;
% QPSK --------------------------------------------------------------------
wifi_params.modulation(2).M = 2; % QPSK
wifi_params.modulation(2).k = 2.^wifi_params.modulation(2).M; % number of constellation points
wifi_params.modulation(2).customSymbolMapping = [3 1 0 2];
wifi_params.modulation(2).hMod = comm.PSKModulator(...
    'ModulationOrder',wifi_params.modulation(2).k,...
    'PhaseOffset',pi/4,...
    'BitInput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(2).customSymbolMapping);
wifi_params.modulation(2).hDemod = comm.PSKDemodulator(...
    'ModulationOrder',wifi_params.modulation(2).k,...
    'PhaseOffset',pi/4,...
    'BitOutput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(2).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(2).hDemodSSD = comm.SphereDecoder(...
    'Constellation',constellation(wifi_params.modulation(2).hDemod),...
    'BitTable',de2bi(wifi_params.modulation(2).hDemod.CustomSymbolMapping,log2(wifi_params.modulation(2).hDemod.ModulationOrder),'left-msb'),...
    'DecisionType','Soft');
wifi_params.modulation(2).hMod_factor = 1;
% 8PSK --------------------------------------------------------------------
wifi_params.modulation(3).M = 3; % 8PSK
wifi_params.modulation(3).k = 2.^wifi_params.modulation(3).M; % number of constellation points
wifi_params.modulation(3).customSymbolMapping = [0 1 3 2 6 7 5 4];
wifi_params.modulation(3).hMod = comm.PSKModulator(...
    'ModulationOrder',wifi_params.modulation(3).k,...
    'PhaseOffset',0,...
    'BitInput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',...
    wifi_params.modulation(3).customSymbolMapping);
wifi_params.modulation(3).hDemod = comm.PSKDemodulator(...
    'ModulationOrder',wifi_params.modulation(3).k,...
    'PhaseOffset',0,...
    'BitOutput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',...
    wifi_params.modulation(3).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(3).hDemodSSD = comm.SphereDecoder(...
    'Constellation',constellation(wifi_params.modulation(3).hDemod),...
    'BitTable',de2bi(wifi_params.modulation(3).hDemod.CustomSymbolMapping,log2(wifi_params.modulation(3).hDemod.ModulationOrder),'left-msb'),...
    'DecisionType','Soft');
wifi_params.modulation(3).hMod_factor = 1;
% 16QAM --------------------------------------------------------------------
wifi_params.modulation(4).M = 4; % 16QAM
wifi_params.modulation(4).k = 2.^wifi_params.modulation(4).M; % number of constellation points
wifi_params.modulation(4).customSymbolMapping = [2 3 1 0 6 7 5 4 ...
    14 15 13 12 10 11 9 8];
wifi_params.modulation(4).hMod = comm.RectangularQAMModulator(...
    'ModulationOrder',wifi_params.modulation(4).k,...
    'NormalizationMethod','Average power',...
    'BitInput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(4).customSymbolMapping);
wifi_params.modulation(4).hDemod = comm.RectangularQAMDemodulator(...
    'BitOutput',true,...
    'ModulationOrder',wifi_params.modulation(4).k,...
    'NormalizationMethod','Average power',...
    'DecisionMethod','Log-likelihood ratio',...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',wifi_params.modulation(4).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(4).hDemodSSD = comm.SphereDecoder(...
    'Constellation',constellation(wifi_params.modulation(4).hDemod),...
    'BitTable',de2bi(wifi_params.modulation(4).hDemod.CustomSymbolMapping,log2(wifi_params.modulation(4).hDemod.ModulationOrder),'left-msb'),...
    'DecisionType','Soft');
wifi_params.modulation(4).hMod_factor = 1/sqrt(10);
% 64QAM --------------------------------------------------------------------
wifi_params.modulation(6).M = 6; % 64QAM
wifi_params.modulation(6).k = 2.^wifi_params.modulation(6).M; % number of constellation points
wifi_params.modulation(6).customSymbolMapping = [4 5 7 6 2 3 1 0 ...
    12 13 15 14 10 11 9 8 ...
    28 29 31 30 26 27 25 24 ...
    20 21 23 22 18 19 17 16 ...
    52 53 55 54 50 51 49 48 ...
    60 61 63 62 58 59 57 56 ...
    44 45 47 46 42 43 41 40 ...
    36 37 39 38 34 35 33 32];
wifi_params.modulation(6).hMod = comm.RectangularQAMModulator(...
    'ModulationOrder',wifi_params.modulation(6).k,...
    'NormalizationMethod','Average power',...
    'BitInput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',...
    wifi_params.modulation(6).customSymbolMapping);
wifi_params.modulation(6).hDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder',wifi_params.modulation(6).k,...
    'NormalizationMethod','Average power',...
    'BitOutput',true,...
    'SymbolMapping','Custom',...
    'CustomSymbolMapping',...
    wifi_params.modulation(6).customSymbolMapping,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'VarianceSource','Input port');
wifi_params.modulation(6).hDemodSSD = comm.SphereDecoder(...
    'Constellation',constellation(wifi_params.modulation(6).hDemod),...
    'BitTable',de2bi(wifi_params.modulation(6).hDemod.CustomSymbolMapping,log2(wifi_params.modulation(6).hDemod.ModulationOrder),'left-msb'),...
    'DecisionType','Soft');
wifi_params.modulation(6).hMod_factor = 1/sqrt(42);
% % 256QAM -- for 802dot11ac and ax only -- not defined correctly yet --------------------------------------------------------------------
% wifi_params.modulation(8).M = 8; %
% wifi_params.modulation(8).k = 2.^wifi_params.modulation(8).M; % number of constellation points
% wifi_params.modulation(8).customSymbolMapping = [];
% wifi_params.modulation(8).hMod = comm.RectangularQAMModulator(...
%     'ModulationOrder',wifi_params.modulation(8).k,...
%     'NormalizationMethod','Average power',...
%     'BitInput',true);
% wifi_params.modulation(8).hDemod = comm.RectangularQAMDemodulator(...
%     'ModulationOrder',wifi_params.modulation(8).k,...
%     'NormalizationMethod','Average power',...
%     'BitOutput',true,...
%     'DecisionMethod','Approximate log-likelihood ratio',...
%     'VarianceSource','Input port');
% wifi_params.modulation(8).hMod_factor = 1/sqrt(170);
% % 1024QAM -- for 802dot11ax only -- not defined correctly yet --------------------------------------------------------------------
% wifi_params.modulation(10).M = 10; %
% wifi_params.modulation(10).k = 2.^wifi_params.modulation(10).M; % number of constellation points
% wifi_params.modulation(10).customSymbolMapping = [];
% wifi_params.modulation(10).hMod = comm.RectangularQAMModulator(...
%     'ModulationOrder',wifi_params.modulation(10).k,...
%     'NormalizationMethod','Average power',...
%     'BitInput',true);
% wifi_params.modulation(10).hDemod = comm.RectangularQAMDemodulator(...
%     'ModulationOrder',wifi_params.modulation(10).k,...
%     'NormalizationMethod','Average power',...
%     'BitOutput',true,...
%     'DecisionMethod','Approximate log-likelihood ratio',...
%     'VarianceSource','Input port');
% wifi_params.modulation(10).hMod_factor = 1/sqrt(682); % 1/sqrt(scf), where scf = 2/3*(M-1) for M-QAM modulations

%% 64Non-Uniform Constellation (64-NUC) --------------------------------------------------------------------
% IEEE 802.11ay SC PHY - special case
if wifi_params.general.useNUC == 1 % constellations from TGay document: 11-15/0601r0, Non-uniform constellations for 64QAM, Daniel Schneider (Sony)
    % see: https://mentor.ieee.org/802.11/dcn/15/11-15-0601-00-00ay-non-uniform-constellations-for-64qam.pptx
    % NUC for CR = 7/8 is not available
    if wifi_params.MCS(i_mcs).CR == 7/8
        warning('NUC for CR = 7/8 is not exactly according to draft');
    end
    % ---------------------------------------------------------------------
    NUC(1).CR = 1/2;
    NUC(1).CR_string = '1/2';
    importNUCcr12;
    NUC(1).SymbolAlphabet = nonUniformConst64cr12(:,3);
    NUC(1).bittable = de2bi(nonUniformConst64cr12(:,1), 6, 'left-msb').';
    clear nonUniformConst64cr12;
    NUC(1).hNUCMod = comm.GeneralQAMModulator(...
        'Constellation',NUC(1).SymbolAlphabet);
    NUC(1).hNUCDemod = comm.GeneralQAMDemodulator(...
        'Constellation',NUC(1).SymbolAlphabet,...
        'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'VarianceSource','Input port');
    % ---------------------------------------------------------------------
    NUC(2).CR = 5/8;
    NUC(2).CR_string = '5/8';
    importNUCcr58;
    NUC(2).SymbolAlphabet = nonUniformConst64cr58(:,3);
    NUC(2).bittable = de2bi(nonUniformConst64cr58(:,1), 6, 'left-msb').';
    clear nonUniformConst64cr58;
    NUC(2).hNUCMod = comm.GeneralQAMModulator(...
        'Constellation',NUC(2).SymbolAlphabet);
    NUC(2).hNUCDemod = comm.GeneralQAMDemodulator(...
        'Constellation',NUC(2).SymbolAlphabet,...
        'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'VarianceSource','Input port');
    % ---------------------------------------------------------------------    
    NUC(3).CR = 3/4;
    NUC(3).CR_string = '3/4';
    importNUCcr34;
    NUC(3).SymbolAlphabet = nonUniformConst64cr34(:,3);
    NUC(3).bittable = de2bi(nonUniformConst64cr34(:,1), 6, 'left-msb').';
    clear nonUniformConst64cr34;
    NUC(3).hNUCMod = comm.GeneralQAMModulator(...
        'Constellation',NUC(3).SymbolAlphabet);
    NUC(3).hNUCDemod = comm.GeneralQAMDemodulator(...
        'Constellation',NUC(3).SymbolAlphabet,...
        'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'VarianceSource','Input port');
    % ---------------------------------------------------------------------
    NUC(4).CR = 13/16;
    NUC(4).CR_string = '13/16';
    importNUCcr1316;
    NUC(4).SymbolAlphabet = nonUniformConst64cr1316(:,3);
    NUC(4).bittable = de2bi(nonUniformConst64cr1316(:,1), 6, 'left-msb').';
    clear nonUniformConst64cr1316;
    NUC(4).hNUCMod = comm.GeneralQAMModulator(...
        'Constellation',NUC(4).SymbolAlphabet);
    NUC(4).hNUCDemod = comm.GeneralQAMDemodulator(...
        'Constellation',NUC(4).SymbolAlphabet,...
        'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'VarianceSource','Input port');
    % ---------------------------------------------------------------------
    NUC(5).CR = 7/8;
    NUC(5).CR_string = '7/8';
    importNUCcr1316;
    NUC(5).SymbolAlphabet = nonUniformConst64cr1316(:,3);
    NUC(5).bittable = de2bi(nonUniformConst64cr1316(:,1), 6, 'left-msb').';
    clear nonUniformConst64cr1316;
    NUC(5).hNUCMod = comm.GeneralQAMModulator(...
        'Constellation',NUC(4).SymbolAlphabet);
    NUC(5).hNUCDemod = comm.GeneralQAMDemodulator(...
        'Constellation',NUC(4).SymbolAlphabet,...
        'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'VarianceSource','Input port');
    % ---------------------------------------------------------------------
    %     nuc_const_all = [cr12, cr58, cr34, cr1316];
    %     nuc_title_all = {'1/2', '5/8', '3/4', '13/16'};
    %     figure;
    %     for inuc = 1:4
    %         subplot(1,4,inuc)
    %         plot(nuc_const_all(:,inuc),'.');
    %         axis equal, axis square
    %         grid on
    %         xlabel('Re'), ylabel('Im')
    %         title(['CR = ', nuc_title_all{inuc}]);
    %     end
    %     print('nuc64','-depsc')
    %     system(['epstopdf nuc64.eps'])
    wifi_params.NUC = NUC;
end

%% Channel coding parameters
m_stbc = 1; % no STBC is used
wifi_params.coding = [];
switch wifi_params.MCS(i_mcs).coding_type
    case 'BCC'
        wifi_params.coding = BCC_params(wifi_params, m_stbc, i_mcs);
    case 'LDPC'
        wifi_params.coding = LDPC_params(wifi_params, m_stbc, i_mcs);
    case {'RS+Blck', 'RS+SPC'}
        wifi_params.coding = RS_Blck_SPC_params(wifi_params, m_stbc, i_mcs);
    case 'none'
        % do nothing (test case)
    otherwise
        error('Unsupported type of channel coding.');
end

if ~(strcmp(wifi_params.general.standard,'802dot11ad') || strcmp(wifi_params.general.standard,'802dot11ay'))
    wifi_params.mapping.n_data_symbols_per_frame = wifi_params.coding.N_sym; % N_SYM
    %     wifi_params.mapping.n_data_symbols_per_frame = 1; % N_SYM
else
    disp('TBD: at load_wifi_params.m, line 442');
end
wifi_params.coding.decision_type = decision_type;

%% Scrambling parameters
wifi_params.scrambling.scr_seed = randi([0 1],1,7); % 7-bit of length, seed selected in a pseudorandom fashion
% AD: see std. IEEE 802.11-2016, ch. 20.3.9, page 2451

%% Interleaving parameters
if strcmp(wifi_params.general.standard,'802dot11ay')
    wifi_params.interleaving = [];
else
    wifi_params.interleaving = Interleaving_params(wifi_params.MCS(i_mcs).coding_type,wifi_params,i_mcs);
end
%% Spreading (IEEE 802.11ad)
wifi_params.spreading.Golay_Seq = Spreading_params;

%% Mapping
wifi_params.general.PHYlayer = wifi_params.MCS(i_mcs).phy_type;

switch wifi_params.general.PHYlayer
    case 'OFDM'
        % it will be necessary in case of change pilot subcarriers position
        map_tmp = false(wifi_params.mapping.n_fft,1);
        % logical map of data
        wifi_params.mapping.map_data = map_tmp;
        wifi_params.mapping.map_data(wifi_params.mapping.ir_data) = true;
        % logical map of pilots
        wifi_params.mapping.map_pilots = map_tmp;
        wifi_params.mapping.map_pilots(wifi_params.mapping.ir_pilots) = true;
        % logical map of DC
        wifi_params.mapping.map_DC = map_tmp;
        wifi_params.mapping.map_DC(wifi_params.mapping.ir_dc) = true;
        % logical phase shift mapping
        wifi_params.mapping.map_phase_rotation = map_tmp;
        wifi_params.mapping.map_phase_rotation(wifi_params.mapping.ir_phase_rotated_sc) = true;
        
        clear map_tmp;
    case 'SC'
        
    case 'LPSC'
        
    case 'Ctrl'
        
    case 'SCraw'
        %         test case, do nothing
    otherwise
        error('Wrong physical layer type (unsupported by simulator)');
end

if strcmp(wifi_params.general.standard,'802dot11ad') || strcmp(wifi_params.general.standard,'802dot11ay') % TBD
    wifi_params.Fs = 2640e6;
else
    wifi_params.Fs = wifi_params.mapping.n_fft*wifi_params.mapping.df;
end

%% GUARD INTERVAL
switch wifi_standard %=====================================================
    case '802dot11g' %-----------------------------------------------------
        switch GuardInterval
            case 'normal'
                wifi_params.cprefix.t_cp = 800e-9;
            otherwise
                error('In IEEE 802.11g ''normal'' cyclic prefix length (800 ns) is valid only');
        end
    case {'802dot11n', '802dot11ac'} %-------------------------------------
        switch GuardInterval
            case 'short'
                wifi_params.cprefix.t_cp = 400e-9;
            case 'normal'
                wifi_params.cprefix.t_cp = 800e-9;
            otherwise
                error('In IEEE 802.11n/ac ''normal'' or ''short'' cyclic prefix length (800 ns or 400 ns) are valid only');
        end
    case '802dot11ax'%-----------------------------------------------------
        switch GuardInterval
            case 'normal'
                wifi_params.cprefix.t_cp = 800e-9;
            case 'doubled'
                wifi_params.cprefix.t_cp = 1600e-9;
            case 'enhanced'
                wifi_params.cprefix.t_cp = 3200e-9;
            otherwise
                error('In IEEE 802.11ax ''normal'', ''doubled'' or ''enhanced'' cyclic prefix length (800 ns, 1600 ns or 3200 ns) are valid only');
        end
    case '802dot11ah' %----------------------------------------------------
        switch GuardInterval
            case 'normal'
                wifi_params.cprefix.t_cp = 8e-6;
            case 'short'
                wifi_params.cprefix.t_cp = 4e-6;
            otherwise
                error('In IEEE 802.11ah ''normal'' or ''short'' cyclic prefix length (8 us or 4 us) are valid only');
        end
    case '802dot11ad'
        wifi_params.cprefix.t_cp = 48.4e-9;
    case '802dot11ay'
        wifi_params.cprefix.t_cp = 48.4e-9;
    otherwise
        error('Wrong or undefined standard');
end

if strcmp(wifi_params.general.PHYlayer,'OFDM')
    wifi_params.cprefix.n_cp = ceil(wifi_params.cprefix.t_cp*wifi_params.Fs);
    wifi_params.cprefix.indices = wifi_params.mapping.n_fft-wifi_params.cprefix.n_cp+1:wifi_params.mapping.n_fft;
end
%% Framing

% disp(mat2str(wifi_params.cprefix.indices))
% disp('TBD: only for AD, see line: 597')
switch wifi_params.general.PHYlayer
    case {'Ctrl'}
        % Framing
        wifi_params.framing.Names = {...
            'STF';...
            'CEF';...
            'Header+Data';...
            'Beamforming Training'}; % AGC ?
        HeadAndDataNumChips = (wifi_params.coding.L_dpfcw+(wifi_params.coding.L_dpcw*(wifi_params.coding.N_cw-2))+wifi_params.coding.L_dplcw+(wifi_params.coding.N_cw*wifi_params.coding.L_cwd))*length(wifi_params.spreading.Golay_Seq.Ga_32);
        wifi_params.framing.ChipLengths = [50*length(wifi_params.spreading.Golay_Seq.Gb_128); 1152; HeadAndDataNumChips; 0]; % in case of Header, see 20.6.3.1.4 Header encoding and modulation
        ChipMapTmp = false(1, sum(wifi_params.framing.ChipLengths));
        % STF
        wifi_params.framing.ChipMap{1} = ChipMapTmp;
        wifi_params.framing.ChipMap{1}(1:wifi_params.framing.ChipLengths(1)) = true;
        % CEF
        wifi_params.framing.ChipMap{2} = ChipMapTmp;
        wifi_params.framing.ChipMap{2}(wifi_params.framing.ChipLengths(1)+(1:wifi_params.framing.ChipLengths(2))) = true;
        % Header+Data
        
        wifi_params.framing.ChipMap{3} = ChipMapTmp;
        wifi_params.framing.ChipMap{3}(wifi_params.framing.ChipLengths(1)+wifi_params.framing.ChipLengths(2)+(1:wifi_params.framing.ChipLengths(3))) = true;
        % Beamforming training
        wifi_params.framing.ChipMap{4} = ChipMapTmp;
        clear ChipMapTmp HeadAndDataNumChips;
        
        wifi_params.framing.Header.Items = {...
            'Differential Encoder Initialization';...
            'Scrambler Initialization';...
            'LENGTH';...
            'Packet Type';...
            'Training Length';...
            'Turnaround';...
            'Reserved';...
            'HCS'};
        wifi_params.framing.Header.ItemsBitLengths = [1, 4, 10, 1, 5, 1, 2, 16];
        wifi_params.framing.Header.ItemsBitLengthsAll = sum(wifi_params.framing.Header.ItemsBitLengths);
        wifi_params.framing.Header.ItemsOrder = 'LSB-first';
        wifi_params.framing.Header.coding = wifi_params.coding; %LDPC_params(wifi_params, 1, 5);
        wifi_params.framing.Header.hCRCgen = comm.CRCGenerator('Polynomial',[16 12 5 0]); %'InitialConditions',ones(1,16)); % not the same as in 15.3.3.7, page 2230
        
        if wifi_params.general.sendDataOnly == 1
            error('Not defined yet for Ctrl PHY');
        end
        
    case {'SC'}
        % Framing
        wifi_params.framing.Names = {...
            'STF';...
            'CEF';...
            'Header*';...
            'Data*';...
            'Beamforming Training'};
        wifi_params.framing.ChipLengths = [...
            2176;...
            1152;...
            (2*wifi_params.mapping.n_tot)+(2*wifi_params.mapping.N_GI);...
            ((wifi_params.coding.N_blks*448/N_ss)+((wifi_params.coding.N_blks/N_ss+1)*wifi_params.mapping.N_GI));...
            0]; % in case of Header, see 20.6.3.1.4 Header encoding and modulation
        ChipMapTmp = false(1, sum(wifi_params.framing.ChipLengths));
        % STF
        wifi_params.framing.ChipMap{1} = ChipMapTmp;
        wifi_params.framing.ChipMap{1}(1:wifi_params.framing.ChipLengths(1)) = true;
        % CEF
        wifi_params.framing.ChipMap{2} = ChipMapTmp;
        wifi_params.framing.ChipMap{2}(wifi_params.framing.ChipLengths(1)+(1:wifi_params.framing.ChipLengths(2))) = true;
        % Header
        wifi_params.framing.ChipMap{3} = ChipMapTmp;
        wifi_params.framing.ChipMap{3}(sum(wifi_params.framing.ChipLengths(1:2))+(1:wifi_params.framing.ChipLengths(3))) = true;
        % Data
        wifi_params.framing.ChipMap{4} = ChipMapTmp;
        wifi_params.framing.ChipMap{4}(sum(wifi_params.framing.ChipLengths(1:3))+(1:wifi_params.framing.ChipLengths(4))) = true;
        % Beamforming training
        wifi_params.framing.ChipMap{5} = ChipMapTmp; % none
        
        wifi_params.framing.Header.Items = {...
            'Scrambler Initialization';...
            'MCS';...
            'LENGTH';...
            'Additional PPDU';...
            'Packet Type';...
            'Training Length';...
            'Aggregation';...
            'Beam Tracking Request';...
            'Last RSSI';...
            'SIFS Response';...
            'Reserved';...
            'HCS'};
        wifi_params.framing.Header.ItemsBitLengths = [7, 5, 18, 1, 1, 5, 1, 1, 4, 1, 4, 16];
        wifi_params.framing.Header.ItemsBitLengthsAll = sum(wifi_params.framing.Header.ItemsBitLengths);
        wifi_params.framing.Header.ItemsOrder = 'LSB-first';
        wifi_params.framing.Header.coding = LDPC_paramsHeader(wifi_params, 1, 5); % LDPC encoder with code rate 3/4
        wifi_params.framing.Header.hCRCgen = comm.CRCGenerator('Polynomial',[16 12 5 0]); %'InitialConditions',ones(1,16)); % not the same as in 15.3.3.7, page 2230
        
        if wifi_params.general.sendDataOnly == 1
            wifi_params.framing.Names = {...
                'Data*'}; % AGC ?
            wifi_params.framing.ChipLengths = [...
                ((wifi_params.coding.N_blks*448)+((wifi_params.coding.N_blks+1)*wifi_params.mapping.N_GI))]; %
            ChipMapTmp = true(1, sum(wifi_params.framing.ChipLengths));
            % Data
            wifi_params.framing.ChipMap = ChipMapTmp;
            clear ChipMapTmp HeadAndDataNumChips;
        end
        
    case {'LPSC'}
        % Framing
        wifi_params.framing.Names = {...
            'STF';...
            'CEF';...
            'Header*';...
            'Data*';...
            'Beamforming Training'};
        wifi_params.framing.ChipLengths = [...
            2176;...
            1152;...
            (2*wifi_params.mapping.n_tot)+(2*wifi_params.mapping.N_GI);...
            (wifi_params.coding.Blck.N_blks*512)+wifi_params.mapping.N_GI;...
            0]; % no Beamforming is used, yet
        ChipMapTmp = false(1, sum(wifi_params.framing.ChipLengths));
        % STF
        wifi_params.framing.ChipMap{1} = ChipMapTmp;
        wifi_params.framing.ChipMap{1}(1:wifi_params.framing.ChipLengths(1)) = true;
        % CEF
        wifi_params.framing.ChipMap{2} = ChipMapTmp;
        wifi_params.framing.ChipMap{2}(wifi_params.framing.ChipLengths(1)+(1:wifi_params.framing.ChipLengths(2))) = true;
        % Header
        wifi_params.framing.ChipMap{3} = ChipMapTmp;
        wifi_params.framing.ChipMap{3}(sum(wifi_params.framing.ChipLengths(1:2))+(1:wifi_params.framing.ChipLengths(3))) = true;
        % Data
        wifi_params.framing.ChipMap{4} = ChipMapTmp;
        wifi_params.framing.ChipMap{4}(sum(wifi_params.framing.ChipLengths(1:3))+(1:wifi_params.framing.ChipLengths(4))) = true;
        % Beamforming training
        wifi_params.framing.ChipMap{5} = ChipMapTmp; % none
        
        wifi_params.framing.Header.Items = {...
            'Scrambler Initialization';...
            'MCS';...
            'LENGTH';...
            'Additional PPDU';...
            'Packet Type';...
            'Training Length';...
            'Aggregation';...
            'Beam Tracking Request';...
            'Last RSSI';...
            'SIFS Response';...
            'Reserved';...
            'HCS'};
        wifi_params.framing.Header.ItemsBitLengths = [7, 5, 18, 1, 1, 5, 1, 1, 4, 1, 4, 16];
        wifi_params.framing.Header.ItemsBitLengthsAll = sum(wifi_params.framing.Header.ItemsBitLengths);
        wifi_params.framing.Header.ItemsOrder = 'LSB-first';
        wifi_params.framing.Header.coding = LDPC_params(wifi_params, 1, 5); % LDPC encoder with code rate 3/4
        wifi_params.framing.Header.hCRCgen = comm.CRCGenerator('Polynomial',[16 12 5 0]); %'InitialConditions',ones(1,16)); % not the same as in 15.3.3.7, page 2230
    case 'OFDM'
        % Windowing
        wifi_params.framing.Windowing = []; % TBD
        % Framing
        wifi_params.framing.Names = {...
            'STF';...
            'CEF';...
            'Header*';...
            'Data*';...
            'Beamforming Training'};
        wifi_params.framing.ChipLengths = [...
            2176;...
            1152;...
            wifi_params.mapping.n_fft+wifi_params.cprefix.n_cp;...
            wifi_params.coding.N_sym*(wifi_params.mapping.n_fft+wifi_params.cprefix.n_cp);...
            0]; % in case of Header, see 20.6.3.1.4 Header encoding and modulation
        ChipMapTmp = false(1, sum(wifi_params.framing.ChipLengths));
        % STF
        wifi_params.framing.ChipMap{1} = ChipMapTmp;
        wifi_params.framing.ChipMap{1}(1:wifi_params.framing.ChipLengths(1)) = true;
        % CEF
        wifi_params.framing.ChipMap{2} = ChipMapTmp;
        wifi_params.framing.ChipMap{2}(wifi_params.framing.ChipLengths(1)+(1:wifi_params.framing.ChipLengths(2))) = true;
        % Header
        wifi_params.framing.ChipMap{3} = ChipMapTmp;
        wifi_params.framing.ChipMap{3}(sum(wifi_params.framing.ChipLengths(1:2))+(1:wifi_params.framing.ChipLengths(3))) = true;
        % Data
        wifi_params.framing.ChipMap{4} = ChipMapTmp;
        wifi_params.framing.ChipMap{4}(sum(wifi_params.framing.ChipLengths(1:3))+(1:wifi_params.framing.ChipLengths(4))) = true;
        % Beamforming training
        wifi_params.framing.ChipMap{5} = ChipMapTmp; % none
        
        wifi_params.framing.Header.Items = {...
            'Scrambler Initialization';...
            'MCS';...
            'LENGTH';...
            'Additional PPDU';...
            'Packet Type';...
            'Training Length';...
            'Aggregation';...
            'Beam Tracking Request';...
            'Tone Pairing Type';...
            'DTP Indicator';...
            'Last RSSI';...
            'Turnaround';...
            'Reserved';...
            'HCS'};
        wifi_params.framing.Header.ItemsBitLengths = [7, 5, 18, 1, 1, 5, 1, 1, 1, 1, 4, 1, 2, 16];
        wifi_params.framing.Header.ItemsBitLengthsAll = sum(wifi_params.framing.Header.ItemsBitLengths);
        wifi_params.framing.Header.ItemsOrder = 'LSB-first';
        wifi_params.framing.Header.coding = LDPC_params(wifi_params, 1, 5); % LDPC encoder with code rate 3/4
        wifi_params.framing.Header.hCRCgen = comm.CRCGenerator('Polynomial',[16 12 5 0]); %'InitialConditions',ones(1,16)); % not the same as in 15.3.3.7, page 2230
        
        if wifi_params.general.sendDataOnly == 1
            error('Not defined yet for OFDM PHY');
        end
        
    otherwise
        error('Wrong type of physical layer')
end

%% Show data rates according to selected MCS
if strcmp(wifi_params.general.standard,'802dot11ad') || strcmp(wifi_params.general.standard,'802dot11ay')
    user_data_rate_str = sprintf(' >>> User data rate: %2.2f Mbps ', maxPHYThroughput/1e6);
    wifi_params.data_rate = maxPHYThroughput;
else
    wifi_params.data_rate = (wifi_params.MCS(i_mcs).N_ss*wifi_params.MCS(i_mcs).CR*wifi_params.mapping.n_data*log2(wifi_params.MCS(i_mcs).M))/((wifi_params.mapping.n_fft/wifi_params.Fs)+(wifi_params.cprefix.t_cp));
    wifi_params.data_rate_PSDU = ((LENGTH*8))/(wifi_params.coding.N_sym*((wifi_params.mapping.n_fft/wifi_params.Fs)+(wifi_params.cprefix.t_cp)));
    
    user_data_rate_str = sprintf(' >>> User data rate: %2.2f Mbps ', wifi_params.data_rate/1e6);
    PSDU_data_rate_str = sprintf(' >>> PSDU data rate: %2.2f Mbps ', wifi_params.data_rate_PSDU/1e6);
end

%% other useful definitions
error_val = zeros(N_frames,length(SNR),length(MCSvec)); % initialize error values (rewrite in future)