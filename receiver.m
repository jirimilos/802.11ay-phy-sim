classdef receiver
    % receiver class for IEEE 802.11 tx
    %
    %
    % Author:	Jiri Milos, DREL FEEC BUT, 2018
    %
    
    properties
        nRxAnt = [];         % number receiving antennas
        fc = [];
        scheme = [];        % MCS
        output = [];
        genie = [];
    end
    
    methods
        function obj = receiver(nReceivers,nRxAnt,fc,scheme) % class constructor
            if nargin ~= 0
                obj(nReceivers) = receiver;
                for ii = 1:nReceivers
                    obj(ii).nRxAnt = nRxAnt;    % number of receiving antennas may be different in every receiver
                    obj(ii).fc = fc;            % carrier frequency should be the same for all receivers
                    obj(ii).scheme = scheme;    % MCS should be the same for all transmitter and receivers
                    obj(ii).output = [];
                end
            end
        end
    end
end