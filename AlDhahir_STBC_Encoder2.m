function [TXA_data, TXB_data] = AlDhahir_STBC_Encoder2(NF, SyD, SyPCP, cplx_data)
% date: 23.1.2020
% STBC encoding across TWO data bursts (Al-Dhahir) suitable for SC-FDE.
% "NF" No. of frames (ROWS of "data FRAMES" input matrix )
% "SyPF" No. of symbols per frame (Columns ...)
% "SyPCP"  no. of symbols per Cyclic Prefix
% "cplx_data" 1xNF*SyPF vector of input cplx symbols
% "mode" declares the number of TX antennas. (mode == 1 for 2 Antennas)
% (mode == 2 for 4 antennas usign CIOD STBC scheme according to ... )
% This function turns input CPLX symbol stream into several output cplx
% streams.
% "CP_nGI" implements Cycle prefix if = 1, otherwise it implements Guard
% Interval (Zeros)

% CHANGES:
% #1
% date: 12.5.2020
% author: Jiri Milos
% Guard interval is added outside this function!
% Other STBC modes were removed.


% % Serial to parallel (Put data into FRAMES):
% dataFRAMES = zeros(NF,SyD);
% for i = 1:NF
%     dataFRAMES(i,:) = cplx_data((i-1)*SyD+1:i*SyD);
% end
dataFRAMES = cplx_data;

%   | s0  s1  |
%   |-s1* s0* |
% (Alamouti scheme)

TXA_data = zeros(NF,SyD); % TX Antenna A
TXB_data = zeros(NF,SyD); % TX Antenna B

for i = 1:2:NF-1
    s0 = dataFRAMES(i,:);
    s1 = dataFRAMES(i+1,:);
    
    aux0 = conj(fliplr(s0)); % s0* time-reversed
    aux1 = conj(fliplr(s1)); % s1* time-reversed
    
    TXA_data(i,:) = s0; % s0(n)
    TXB_data(i,:) = s1; % s1(n)
    
    TXA_data(i+1,:) = -[aux1(end), aux1(1:end-1)];% -s1*((-n)N)
    TXB_data(i+1,:) = [aux0(end), aux0(1:end-1)];% s0*((-n)N)   
end

end

