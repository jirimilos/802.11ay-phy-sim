function output = ldpcEncoder(input,wifi_params, MCS)
%LPDCEncoder Function performs channel encoding using LDPC block coding 
%   Author: Jiri Milos, DREL FEEC BUT, 2020

wifiStandard = wifi_params.general.standard;
phyType = MCS.phy_type;

switch wifiStandard
    case '802dot11ad' % ===================================================
   
        
        
        
        
        
        
    case '802dot11ay' % ===================================================
        
        
        
        
        
        
        
    otherwise % ===========================================================
            error('Wrong IEEE 802.11 standard');
end









output = input;
end

