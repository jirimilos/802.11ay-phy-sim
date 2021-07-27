function PSDU_tx = generatePSDU(wifi_params)

% Number of PSDU data octets
LENGTH = wifi_params.general.LENGTH;
if wifi_params.general.sendAllZeros == true
    PSDU_tx = zeros(LENGTH*8,1);
elseif wifi_params.general.sendAllOnes == true
    PSDU_tx = ones(LENGTH*8,1);
else
    PSDU_tx = randi([0 1],LENGTH*8,1);
end


end

