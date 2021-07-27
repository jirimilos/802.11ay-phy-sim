function [c_outSSD, c_outHD] = DCM_demod(inputModSymbols,H_est, ModulationType)
%DCM_MOD Dual carrier modulation (DCM) demodulator function
% It performs inverse operation to functions described in Standard IEEE 802.11ad in OFDM PHY
% layer in Sec. from 20.5.3.2.4.3 to 20.5.3.2.4.5
% Static tone pairing only
%
% Author: Jiri Milos, DREL FEEC BUT, 2020

switch ModulationType
    case 'QPSK'
        M = 4;
        Qmat = (1/sqrt(5))*[1,2;-2,1];
        QmatR = pinv(Qmat);
        bittable = logical([1 1 0 0;...
            1 0 1 0]);
        SymbolAlphabet = ((1/sqrt(2))*((2*bittable(1,:).'-1) +1j*(2*bittable(2,:).'-1))).';
    case '16QAM'
        M = 16;
        Qmat = 1;
        QmatR = 1;
        bittable = logical([...
            0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1;...
            0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1;...
            0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1;...
            0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1]);
        SymbolAlphabet = (1/sqrt(10))*[...
            -3-1j*3,  -3-1j*1,  -3+1j*3,  -3+1j*1,...
            -1-1j*3,  -1-1j*1,  -1+1j*3,  -1+1j*1,...
            3-1j*3,   3-1j*1,   3+1j*3,   3+1j*1,...
            1-1j*3,   1-1j*1,   1+1j*3,   1+1j*1];
    case '64QAM'
        M = 64;
        Qmat = 1;
        QmatR = 1;
        bittable = logical([...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;...
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;...
            0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1;...
            0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1;...
            0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1;...
            0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1]);
        SymbolAlphabet = (1/sqrt(42))*[...
            -7-1j*7,   -7-1j*5,   -7-1j*1,   -7-1j*3,   -7+1j*7,   -7+1j*5,   -7+1j*1,   -7+1j*3,...
            -5-1j*7,   -5-1j*5,   -5-1j*1,   -5-1j*3,   -5+1j*7,   -5+1j*5,   -5+1j*1,   -5+1j*3,...
            -1-1j*7,   -1-1j*5,   -1-1j*1,   -1-1j*3,   -1+1j*7,   -1+1j*5,   -1+1j*1,   -1+1j*3,...
            -3-1j*7,   -3-1j*5,   -3-1j*1,   -3-1j*3,   -3+1j*7,   -3+1j*5,   -3+1j*1,   -3+1j*3,...
            7-1j*7,    7-1j*5,    7-1j*1,    7-1j*3,    7+1j*7,    7+1j*5,    7+1j*1,    7+1j*3,...
            5-1j*7,    5-1j*5,    5-1j*1,    5-1j*3,    5+1j*7,    5+1j*5,    5+1j*1,    5+1j*3,...
            1-1j*7,    1-1j*5,    1-1j*1,    1-1j*3,    1+1j*7,    1+1j*5,    1+1j*1,    1+1j*3,...
            3-1j*7,    3-1j*5,    3-1j*1,    3-1j*3,    3+1j*7,    3+1j*5,    3+1j*1,    3+1j*3];
    otherwise
        error('wrong modulation type');
end

mR = length(inputModSymbols);
if strcmp(ModulationType,'64QAM') == 0
    dR = [inputModSymbols(1); inputModSymbols(2)];
    % Get QPSK or 16QAM modulation symbols
    xR = QmatR*dR;
    
    [Q1,R1] = qr(H_est(1)); 
    [Q2,R2] = qr(H_est(2)); 
    
    
    % symbols to demodulate
    rx_user_sym_tmp1 = xR(1);
    rx_user_sym_tmp2 = xR(2);
    % Zero forcing equalization
    rx_layer_x1 = (H_est(1).^(-1)).*rx_user_sym_tmp1; % ZF % var. 1
    rx_layer_x2 = (H_est(2).^(-1)).*rx_user_sym_tmp2; % ZF % var. 2
    
    % soft demodulation
    c_out1_SSD = LTE_softsphere(rx_layer_x1, rx_user_sym_tmp1, Q1, R1, SymbolAlphabet, (bittable), 1, log2(M));
    c_out2_SSD = LTE_softsphere(rx_layer_x2, rx_user_sym_tmp2, Q2, R2, SymbolAlphabet, (bittable), 1, log2(M));
else
    dR = inputModSymbols(1);
    % Get QPSK or 16QAM modulation symbols
    xR = QmatR*dR;
    
    [Q,R] = qr(H_est); % for AWGN only
    
    % symbol to demodulate
    rx_user_sym_tmp = xR;
    % Zero forcing equalization
    rx_layer_x = (H_est.^(-1)).*rx_user_sym_tmp; % ZF % var. 1
    % soft demodulation
    c_out1_SSD = LTE_softsphere(rx_layer_x, rx_user_sym_tmp, Q, R, SymbolAlphabet, (bittable), 1, log2(M));
    
end

switch ModulationType
    case 'QPSK'
        % arrange output LLR values
        c_outSSD = [c_out1_SSD(1), c_out2_SSD(1), c_out1_SSD(2), c_out2_SSD(2)];
    case '16QAM'
        % arrange output LLR values
        c_outSSD = [c_out1_SSD; c_out2_SSD];
    case '64QAM'
        % arrange output LLR values
        c_outSSD = c_out1_SSD;
    otherwise
        error('Wrong DCM modulation type.')
end

c_outHD = double(c_outSSD>0);
end

