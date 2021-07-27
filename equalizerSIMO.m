function [output, EqNom, EqDenom] = equalizerSIMO(inputPerRxAnt, hDPerRxAnt, noiseVar, EqMode)
%equalizerSIMO Funkce pro ekvalizaci prijatych symbolu (SIMO model)
% Maximum Ratio Combining

hD = squeeze(hDPerRxAnt);
input = squeeze(inputPerRxAnt);

switch EqMode
    case 'ZF'
        EqNom = conj(hD);
        EqDenom = conj(hD).*hD; % Zero Forcing
    case 'MMSE'
        EqNom = conj(hD);
        EqDenom = (conj(hD).*hD)+noiseVar; % MMSE
    otherwise
        error('Only ZF or MMSE equalizer is available');
end

output = sum(input.*EqNom,2)./sum(EqDenom,2);
end