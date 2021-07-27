function [output, Eq] = equalizerSISO(input, hD, noiseVar, EqMode)
%equalizerSISO Funkce pro ekvalizaci prijatych symbolu (SISO model)

switch EqMode
    case 'ZF'
        Eq = (conj(hD))./((conj(hD)).*hD); % Zero Forcing
    case 'MMSE'
        Eq = (conj(hD))./((conj(hD).*hD)+noiseVar); % MMSE
    otherwise
        error('Only ZF or MMSE equalizer is available');
end

output = input.*Eq;
end

