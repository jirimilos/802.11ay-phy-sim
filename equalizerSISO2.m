function [outputBlocks] = equalizerSISO2(inputBlocks, H, noiseVar, SNR_dB)
%equalizerSISO Funkce pro ekvalizaci prijatych symbolu (SISO model)

SNR_lin = 10^(SNR_dB/10); % [-]


N_blks = size(inputBlocks,1);
Hmat = repmat(H, N_blks, 1);

% Y_rx_eq = (fft(inputBlocks.').').*conj(Hmat)./((abs(Hmat).^2)+(1/SNR_lin));
Y_rx_eq = (fft(inputBlocks.').').*conj(Hmat)./((abs(Hmat).^2)+noiseVar);
outputBlocks = ifft(Y_rx_eq.').';

end

