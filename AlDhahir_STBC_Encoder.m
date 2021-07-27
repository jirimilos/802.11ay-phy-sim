function [STBC_out] = AlDhahir_STBC_Encoder(cplx_data, wifi_params)
% function [STBC_data] = AlDhahir_STBC_Encoder(N_blks, SyPF, SyPCP, cplx_data, mode, CP_nGI)
% date: 23.1.2020
% STBC encoding across TWO data bursts (Al-Dhahir) suitable for SC-FDE.
% "N_blks" No. of frames (ROWS of "data FRAMES" input matrix )
% "SyPF" No. of symbols per frame (Columns ...)
% "SyPCP"  no. of symbols per Cyclic Prefix
% "cplx_data" 1xN_blks*SyPF vector of input cplx symbols
% "mode" declares the number of TX antennas. (mode == 1 for 2 Antennas)
% (mode == 2 for 4 antennas usign CIOD STBC scheme according to ... )
% This function turns input CPLX symbol stream into several output cplx
% streams.
% "CP_nGI" implements Cycle prefix if = 1, otherwise it implements Guard
% Interval (Zeros)

% Author: Radim Zedka
% Changes: Jiri Milos, 2020/03/20

mimoMode = wifi_params.antConfig.mode;
nTxAnt = wifi_params.antConfig.nTxAnt;

mode = 1; % test, 2x1 and 2x2 only

GIType = wifi_params.general.guardInterval.type;
Ga64 = wifi_params.spreading.Golay_Seq.Ga_64;
Gb64 = wifi_params.spreading.Golay_Seq.Gb_64;
N_blks = wifi_params.coding.N_blks; % N_blks = N_blks

n_fft = wifi_params.mapping.n_fft;

SyPF = 448; % set directly
% SyPCP = 64;
%
% cplx_data = cplx_data.';
%
% % Serial to parallel (Put data into BLOCKS (=FRAMES)):
% dataFRAMES = zeros(N_blks,SyPF);
% for i = 1:N_blks
%     dataFRAMES(i,:) = cplx_data((i-1)*SyPF+1:i*SyPF);
% end

dataFRAMES = cplx_data;

if (nTxAnt == 2)  % 2 TX antennas
    if strcmp(mimoMode,'TxD') % 2x1 MISO or 2x2 MIMO - Alamouti scheme
        %     (mode == 1)
        %   | s0  s1  |
        %   |-s1* s0* |
        % (Alamouti scheme)
        
%         TXA_data = zeros(N_blks,n_fft); % TX Antenna A
%         TXB_data = zeros(N_blks,n_fft); % TX Antenna B
        for i = 1:2:N_blks-1
            s1 = dataFRAMES(i,:);
            s2 = dataFRAMES(i+1,:);
            
            aux1 = conj(fliplr(s1)); % s0* time-reversed
            aux2 = conj(fliplr(s2)); % s1* time-reversed
            
            TXA_data(i,:) = s1; % s0(n)
            TXB_data(i,:) = s2; % s1(n)
%             % Al-Dhahir
%             TXA_data(i+1,:) = -[aux2(end), aux2(1:end-1)];% -s1*((-n)N)
%             TXB_data(i+1,:) = [aux1(end), aux1(1:end-1)];% s0*((-n)N)
            % UW SC PHY
            TXA_data(i+1,:) = -[aux2];% -s1*((-n-1)N)
            TXB_data(i+1,:) = [aux1];% s0*((-n-1)N)
            
        end
        % Adding Cyclic Prefix and Reshaping back into a cplx vector:
        % Adding the last guard on the end of the final block
        if strcmp(GIType, 'CP') == 1 % Cyclic prefix
            TXA_data_cp = [TXA_data(:,(SyPF-SyPCP+1):SyPF), TXA_data(:,1:SyPF)];
            TXB_data_cp = [TXB_data(:,(SyPF-SyPCP+1):SyPF), TXB_data(:,1:SyPF)]; % ROWs == FRAMEs
            
            TXA_data_append = TXA_data(N_blks,(SyPF-SyPCP+1):SyPF);
            TXB_data_append = TXB_data(N_blks,(SyPF-SyPCP+1):SyPF);
            
            STBC_data = zeros(2, (N_blks*(SyPF + SyPCP))+SyPCP);
            STBC_data(1,:) = [reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXA_data_append];
            STBC_data(2,:) = [reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXB_data_append];
            
        elseif strcmp(GIType, 'GI') == 1 % Guard interval (is included!)
            puvodni = 0;
            if puvodni == 1
                % puvodni verze
                k_Ga64 = 0:length(Ga64)-1;
                Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
                                
                % Ga64 to prepend
                TXA_data_prepend = Ga64_rotated;
                TXB_data_prepend = Ga64_rotated;
                
                TXA_data_cp = [TXA_data_prepend, reshape(TXA_data(:,1:n_fft).',1,[])];
                TXB_data_cp = [TXB_data_prepend, reshape(TXB_data(:,1:n_fft).',1,[])]; % ROWs == FRAMEs
                
                STBC_data = zeros(2, length(TXA_data_cp));
                STBC_data(1,:) = TXA_data_cp;
                STBC_data(2,:) = TXB_data_cp;
            else
                % test verze zohlednujici upravu guard intervalu ve vetvi B (uprava kvuli odstraneni IBI)
                k_Ga64 = 0:length(Ga64)-1;
                Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
                Gb64_rotated = Gb64.*exp(1j*pi*k_Ga64/2);
                
                % TXA and TXB - append Ga64 (always the same)
                TXA_data_Ga64appended = [TXA_data, repmat(Ga64_rotated,N_blks,1)];
                TXB_data_Gb64appended = [TXB_data, repmat(Gb64_rotated,N_blks,1)];
%                 TXB_data_Ga64appended = [TXB_data, repmat(Gb64_rotated,N_blks,1)];
                
                TXA_data_cp = [Ga64_rotated, reshape(TXA_data_Ga64appended(:,1:n_fft).',1,[])];
                TXB_data_cp = [Gb64_rotated, reshape(TXB_data_Gb64appended(:,1:n_fft).',1,[])]; % ROWs == FRAMEs
%                 TXB_data_cp = [Gb64_rotated, reshape(TXB_data_Ga64appended(:,1:n_fft).',1,[])]; % ROWs == FRAMEs
                STBC_data = zeros(2, length(TXA_data_cp));
                STBC_data(1,:) = TXA_data_cp;
                STBC_data(2,:) = TXB_data_cp;
            end
            
        elseif strcmp(GIType, 'zeros')
            TXA_data_cp = [zeros(N_blks, SyPCP), TXA_data(:,1:SyPF)];
            TXB_data_cp = [zeros(N_blks, SyPCP), TXB_data(:,1:SyPF)];
            
            TXA_data_append = zeros(1, SyPCP);
            TXB_data_append = zeros(1, SyPCP);
            
            STBC_data = zeros(2, (N_blks*(SyPF + SyPCP))+SyPCP);
            STBC_data(1,:) = [reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXA_data_append];
            STBC_data(2,:) = [reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXB_data_append];
            
        else
            error('Wrong CP type');
        end
        
        
        STBC_out = STBC_data.';
        
    elseif strcmp(mimoMode,'MIMO') % 2x2 MIMO - Spatial Multiplexing using encoding based on Alamouti scheme (STBC), see: IEEE 802.11-11/1712r00
        
        %   | s0  s2  |
        %   |-s1* s3* |
        % based on (Alamouti scheme)
        
        TXA_data = zeros(N_blks,n_fft); % TX Antenna A
        TXB_data = zeros(N_blks,n_fft); % TX Antenna B
        for i = 1:2:N_blks-1
            s1 = dataFRAMES(i,:);
            s2 = dataFRAMES(i+1,:);
            
            aux1 = conj(fliplr(s1)); % s0* time-reversed
            aux2 = conj(fliplr(s2)); % s1* time-reversed
            
            TXA_data(i,:) = s1; % s0(n)
            TXB_data(i,:) = s2; % s1(n)
            
            TXA_data(i+1,:) = -[aux2(end), aux2(1:end-1)];% -s1*((-n)N)
            TXB_data(i+1,:) = [aux1(end), aux1(1:end-1)];% s0*((-n)N)
            
        end
        % Adding Cyclic Prefix and Reshaping back into a cplx vector:
        % Adding the last guard on the end of the final block
        if strcmp(GIType, 'CP') == 1 % Cyclic prefix
            TXA_data_cp = [TXA_data(:,(SyPF-SyPCP+1):SyPF), TXA_data(:,1:SyPF)];
            TXB_data_cp = [TXB_data(:,(SyPF-SyPCP+1):SyPF), TXB_data(:,1:SyPF)]; % ROWs == FRAMEs
            
            TXA_data_append = TXA_data(N_blks,(SyPF-SyPCP+1):SyPF);
            TXB_data_append = TXB_data(N_blks,(SyPF-SyPCP+1):SyPF);
            
            STBC_data = zeros(2, (N_blks*(SyPF + SyPCP))+SyPCP);
            STBC_data(1,:) = [reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXA_data_append];
            STBC_data(2,:) = [reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXB_data_append];
            
        elseif strcmp(GIType, 'GI') == 1 % Guard interval (is included!)
            k_Ga64 = 0:length(Ga64)-1;
            Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
            
            % Ga64 to prepend
            TXA_data_prepend = Ga64_rotated;
            TXB_data_prepend = Ga64_rotated;
            
            TXA_data_cp = [TXA_data_prepend, reshape(TXA_data(:,1:n_fft).',1,[])];
            TXB_data_cp = [TXB_data_prepend, reshape(TXB_data(:,1:n_fft).',1,[])]; % ROWs == FRAMEs
            
            STBC_data = zeros(2, length(TXA_data_cp));
            STBC_data(1,:) = TXA_data_cp;
            STBC_data(2,:) = TXB_data_cp;
            
        elseif strcmp(GIType, 'zeros')
            TXA_data_cp = [zeros(N_blks, SyPCP), TXA_data(:,1:SyPF)];
            TXB_data_cp = [zeros(N_blks, SyPCP), TXB_data(:,1:SyPF)];
            
            TXA_data_append = zeros(1, SyPCP);
            TXB_data_append = zeros(1, SyPCP);
            
            STBC_data = zeros(2, (N_blks*(SyPF + SyPCP))+SyPCP);
            STBC_data(1,:) = [reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXA_data_append];
            STBC_data(2,:) = [reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP)), TXB_data_append];
            
        else
            error('Wrong CP type');
        end
        
        
        STBC_out = STBC_data.';
        
        
    end
    
    
elseif (nTxAnt == 4)  % CIOD 4x1 (Md. Zafar Ali Khan, B. Sundar Rajan, 2006)
    error()
    %      A   B    C   D
    %   | s0  s1    0   0 |
    %   |-s1* s0*   0   0 |
    %   |  0   0   s2  s3 |
    %   |  0   0  -s3* s2*|
    
    TXA_data = zeros(N_blks, SyPF); % TX Antenna A
    TXB_data = zeros(N_blks, SyPF); % TX Antenna B
    TXC_data = zeros(N_blks, SyPF); % TX Antenna C
    TXD_data = zeros(N_blks, SyPF); % TX Antenna D
    
    for i = 1:4:N_blks-3
        s1 = dataFRAMES(i,:);
        s2 = dataFRAMES(i+1,:);
        s2 = dataFRAMES(i+2,:);
        s3 = dataFRAMES(i+3,:);
        
        aux1 = conj(fliplr(s1));% s0* time-reversed
        aux2 = conj(fliplr(s2));% s1* time-reversed
        aux2 = conj(fliplr(s2));% s2* time-reversed
        aux3 = conj(fliplr(s3));% s3* time-reversed
        
        TXA_data(i,:) = s1; % s0
        TXB_data(i,:) = s2; % s1
        %         TXC_data(i,:) = 0; %
        %         TXD_data(i,:) = 0; %
        
        TXA_data(i+1,:) = -[aux2(end) aux2(1:end-1)]; % -s1*((-n)N)
        TXB_data(i+1,:) =  [aux1(end) aux1(1:end-1)]; % s0*((-n)N)
        %         TXC_data(i+1,:) = 0; %
        %         TXD_data(i+1,:) = 0; %
        
        %         TXA_data(i+2,:) = 0; %
        %         TXB_data(i+2,:) = 0; %
        TXC_data(i+2,:) = s2; % s2
        TXD_data(i+2,:) = s3; % s3
        
        %         TXA_data(i+3,:) = 0; %
        %         TXB_data(i+3,:) = 0; %
        TXC_data(i+3,:) = -[aux3(end) aux3(1:end-1)]; % -s3*
        TXD_data(i+3,:) =  [aux2(end) aux2(1:end-1)]; % s2*
    end
    % Adding Cyclic Prefix and Reshaping back into a cplx vector:
    TXA_data_cp = [TXA_data(:,(SyPF-SyPCP+1):SyPF) TXA_data(:,1:SyPF)];
    TXB_data_cp = [TXB_data(:,(SyPF-SyPCP+1):SyPF) TXB_data(:,1:SyPF)]; % ROWs == FRAMEs
    TXC_data_cp = [TXC_data(:,(SyPF-SyPCP+1):SyPF) TXC_data(:,1:SyPF)];
    TXD_data_cp = [TXD_data(:,(SyPF-SyPCP+1):SyPF) TXD_data(:,1:SyPF)]; % ROWs == FRAMEs
    
    % Reshape into (4 x N_blks*(SyPF+SyPCP)) Matrix:
    STBC_data = zeros(4, N_blks*(SyPF + SyPCP));
    STBC_data(1,:) = reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(2,:) = reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(3,:) = reshape(TXC_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(4,:) = reshape(TXD_data_cp.', 1, N_blks*(SyPF+SyPCP));
    
    
elseif (nTxAnt == 2) && (mode == 3) % 4 TX antennas (4x1 MISO Jafarkhani Quasi-Orthogonal)
    
    %   |  s0   s1   s2  s3  |
    %   | -s1*  s0* -s3* s2* |
    %   | -s2* -s3*  s0* s1* |
    %   |  s3  -s2  -s1  s0  |
    TXA_data = zeros(N_blks, SyPF); % TX Antenna A
    TXB_data = zeros(N_blks, SyPF); % TX Antenna B
    TXC_data = zeros(N_blks, SyPF); % TX Antenna C
    TXD_data = zeros(N_blks, SyPF); % TX Antenna D
    
    for i = 1:4:N_blks-3
        s1 = dataFRAMES(i,:);
        s2 = dataFRAMES(i+1,:);
        s2 = dataFRAMES(i+2,:);
        s3 = dataFRAMES(i+3,:);
        
        
        aux1 = conj(fliplr(s1));% s0* time-reversed
        aux2 = conj(fliplr(s2));% s1* time-reversed
        aux2 = conj(fliplr(s2));% s2* time-reversed
        aux3 = conj(fliplr(s3));% s3* time-reversed
        
        
        TXA_data(i,:) = s1; % s0
        TXB_data(i,:) = s2; % s1
        TXC_data(i,:) = s2; % s2
        TXD_data(i,:) = s3; % s3
        
        TXA_data(i+1,:) = -[aux2(end) aux2(1:end-1)]; % -s1*((-n)N)
        TXB_data(i+1,:) = [aux1(end) aux1(1:end-1)]; % s0*((-n)N)
        TXC_data(i+1,:) = -[aux3(end) aux3(1:end-1)]; % -s3*((-n)N)
        TXD_data(i+1,:) = [aux2(end) aux2(1:end-1)]; % s2*((-n)N)
        
        TXA_data(i+2,:) = -[aux2(end) aux2(1:end-1)]; % -s2*((-n)N)
        TXB_data(i+2,:) = -[aux3(end) aux3(1:end-1)]; % -s3*((-n)N)
        TXC_data(i+2,:) = [aux1(end) aux1(1:end-1)]; % s0*((-n)N)
        TXD_data(i+2,:) = [aux2(end) aux2(1:end-1)]; % s1*((-n)N)
        
        TXA_data(i+3,:) = s3; % s3
        TXB_data(i+3,:) = -s2; % -s2
        TXC_data(i+3,:) = -s2; % -s1
        TXD_data(i+3,:) = s1; % s0
    end
    % Adding Cyclic Prefix and Reshaping back into a cplx vector:
    TXA_data_cp = [TXA_data(:,(SyPF-SyPCP+1):SyPF) TXA_data(:,1:SyPF)];
    TXB_data_cp = [TXB_data(:,(SyPF-SyPCP+1):SyPF) TXB_data(:,1:SyPF)]; % ROWs == FRAMEs
    TXC_data_cp = [TXC_data(:,(SyPF-SyPCP+1):SyPF) TXC_data(:,1:SyPF)];
    TXD_data_cp = [TXD_data(:,(SyPF-SyPCP+1):SyPF) TXD_data(:,1:SyPF)]; % ROWs == FRAMEs
    
    % Reshape into (4 x N_blks*(SyPF+SyPCP)) Matrix:
    STBC_data = zeros(4, N_blks*(SyPF + SyPCP));
    STBC_data(1,:) = reshape(TXA_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(2,:) = reshape(TXB_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(3,:) = reshape(TXC_data_cp.', 1, N_blks*(SyPF+SyPCP));
    STBC_data(4,:) = reshape(TXD_data_cp.', 1, N_blks*(SyPF+SyPCP));
    
    
end

end