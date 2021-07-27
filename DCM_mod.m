function [d_out] = DCM_mod(inputBitWord,ModulationType)
%DCM_MOD Dual carrier modulation (DCM) modulator function
% It performs operation according to Standard IEEE 802.11ad in OFDM PHY
% layer. See Sec. from 20.5.3.2.4.3 to 20.5.3.2.4.5
% Static tone pairing only
%
% Author: Jiri Milos, DREL FEEC BUT, 2020
m = length(inputBitWord);

switch ModulationType
    case 'QPSK'
        Qmat = (1/sqrt(5))*[1,2;-2,1];
        c02 = [inputBitWord(1), inputBitWord(3)];
        c13 = [inputBitWord(2), inputBitWord(4)];
        % Produce QPSK modulation symbols
        x0 = (1/sqrt(2))*((2*c02(1)-1)+1j*(2*c02(2)-1));
        x1 = (1/sqrt(2))*((2*c13(1)-1)+1j*(2*c13(2)-1));
    case '16QAM'
        Qmat = 1;
        c0123 = inputBitWord(1:4);
        c4567 = inputBitWord(5:8);
        % Produce 16QAM modulation symbols
        x0 = (1/sqrt(10))*(((4*c0123(1)-2)-((2*c0123(1)-1).*(2*c0123(2)-1)))+1j*((4*c0123(3)-2)-((2*c0123(3)-1).*(2*c0123(4)-1))));
        x1 = (1/sqrt(10))*(((4*c4567(1)-2)-((2*c4567(1)-1).*(2*c4567(2)-1)))+1j*((4*c4567(3)-2)-((2*c4567(3)-1).*(2*c4567(4)-1))));
    case '64QAM'
        Qmat = 1;
        c012345 = inputBitWord(1:6);
        x0 = (1/sqrt(42))*(((8*c012345(1)-4)-((2*c012345(1)-1).*(4*c012345(2)-2))+((2*c012345(1)-1)).*(2*c012345(2)-1).*(2*c012345(3)-1))+1j*((8*c012345(4)-4)-((2*c012345(4)-1).*(4*c012345(5)-2))+((2*c012345(4)-1)).*(2*c012345(5)-1).*(2*c012345(6)-1)));
        x1 = [];
    otherwise
        error('wrong modulation type');
end

% produce single or two constellation points
x01 = [x0;x1];
% DCM mapping
d = Qmat*x01;
% send to output
d_out = d;
end

