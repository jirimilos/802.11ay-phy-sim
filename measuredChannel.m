classdef measuredChannel
    %measuredChannel class for processing and vizualization of measured traces obtained by indoor measurement using VNA
    % at Department of Radio Electronics, Brno University of Technology
    %
    % Author: Jiri Milos, DREL, 2020
    
    properties
        nPoints;
        fMin;
        fMax;
        fVec;
        nTxAnt;
        nRxAnt;
        Type;
        CIR = [];
        CTF = [];
        bandwidth;
        timeResolution;
        tVec;
        cPower = [];
    end
    
    methods
        function obj = measuredChannel(chanInput,inputType, nPoints, nTxAnt, nRxAnt, fMin, fMax)
            %UNTITLED7 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Type = inputType;
            obj.nPoints = nPoints;
            obj.nTxAnt = nTxAnt;
            obj.nRxAnt = nRxAnt;
            obj.fMin = fMin;
            obj.fMax = fMax;
            obj.fVec = linspace(fMin, fMax, nPoints);
            obj.bandwidth = obj.fMax-obj.fMin;
            obj.timeResolution = 1/obj.bandwidth;
            obj.tVec = 0:obj.timeResolution:obj.timeResolution*(obj.nPoints-1);
            switch obj.Type
                case 'CTF'
                    ctf = zeros(nTxAnt, nPoints, nRxAnt);
                    for i=1:4
                        for j=1:4
                            ctfRealPart = squeeze(real(chanInput(i,j,:)));
                            ctfImagPart = squeeze(imag(chanInput(i,j,:)));
                            ctf(i,:,j) = ctfRealPart+1j*ctfImagPart;
                            %
%                             obj.cPower(i,j).P0 = (1/length(ctf(i,:,j)))*sum(abs(ctf(i,:,j)).^2);
%                             obj.cPower(i,j).P0_rms = sqrt(obj.cPower(i,j).P0);
%                             obj.cPower(i,j).P0_rms_dB = (1/2)*20*log10(obj.cPower(i,j).P0);
                        end
                    end
                    obj.CTF = ctf;
                    obj.CIR = obj.ctf2cir;
                case 'CIR'
                    cir = zeros(nTxAnt, nPoints, nRxAnt);
                    for i=1:nTxAnt
                        for j=1:nRxAnt
                            cirRealPart = squeeze(real(chanInput(i,j,:)));
                            cirImagPart = squeeze(imag(chanInput(i,j,:)));
                            cir(i,:,j) = cirRealPart+1j*cirImagPart;
                        end
                    end
                    obj.CIR = cir;
                    obj.CTF = obj.cir2ctf;
                otherwise
                    error('Nonexistent option');
            end
        end
        
        function hFig = showCir(obj, figureNum, figureType)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(obj.CIR) == 1
                error('Object stores Channel transfer function (CTF) only');
            end
            % create figure
            hFig = figure(figureNum);
            
            yLimits = obj.getYlimits; % obtain yLimits, jeste zkontrolovat
            
            pdp = zeros(size(obj.CIR));
            pdpdB = pdp;
            k = 1;
            for i=1:obj.nTxAnt
                for j=1:obj.nRxAnt
                    pdp(i,:,j) = abs(obj.CIR(i,:,j).*conj(obj.CIR(i,:,j)));
                    pdpdB(i,:,j) = 20*log10(pdp(i,:,j));
                    % switch plot type
                    if strcmp(figureType,'single') == 1
                        plot(obj.tVec, pdpdB(i,:,j));
                        hold all, grid on;
%                         ylim(yLimits(:,:,2))
                    elseif strcmp(figureType,'matrix') == 1
                        subplot(obj.nTxAnt, obj.nRxAnt,k);
                        if i==j
                            plot(obj.tVec, pdpdB(i,:,j),'r');
                        else
                            plot(obj.tVec, pdpdB(i,:,j),'b');
                        end
                        grid on;
                        k = k+1;
                        title(['TX',num2str(i),' ---> RX',num2str(j)]);
%                         ylim(yLimits(:,:,2))
                    else
                        error('Unsupproted option');
                    end
                    xlabel('Time [s]')
                    ylabel('PDP [dB]')
                end
                
            end
            
            
        end
        
        function hFig = showCtf(obj, figureNum, figureType)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(obj.CTF) == 1
                error('Object stores Channel impulse response (CIR) only');
            end
            % create figure
            hFig = figure(figureNum);
            
            yLimits = obj.getYlimits; % obtain yLimits
                        
            ctf = zeros(size(obj.CTF));
            ctfdB = ctf;
            k = 1;
            for i=1:obj.nTxAnt
                for j=1:obj.nRxAnt
                    ctf(i,:,j) = abs(obj.CTF(i,:,j).*conj(obj.CTF(i,:,j)));
                    ctfdB(i,:,j) = 20*log10(ctf(i,:,j));
                    % switch plot type
                    if strcmp(figureType,'single') == 1
                        plot(obj.fVec, ctfdB(i,:,j));
                        hold all, grid on;
%                         ylim(yLimits(:,:,1))
                    elseif strcmp(figureType,'matrix') == 1
                        subplot(obj.nTxAnt, obj.nRxAnt,k);
                        if i==j
                            plot(obj.fVec, ctfdB(i,:,j),'r');
                        else
                            plot(obj.fVec, ctfdB(i,:,j),'b');
                        end
                        grid on;
                        k = k+1;
                        title(['TX',num2str(i),' ---> RX',num2str(j)]);
%                         ylim(yLimits(:,:,1))
                    else
                        error('Unsupproted option');
                    end
                    xlabel('Frequency [Hz]')
                    ylabel('CTF [dB]')
                    
                end
                
            end
            
            
        end
        
        function obj = restrictCtf(obj, fResMin, fResMax)
            if isempty(obj.CTF) == 1
                error('CTF is not defined');
            end
            
            [~, indfResMin] = find(abs(obj.fVec-fResMin) == 0);
            [~, indfResMax] = find(abs(obj.fVec-fResMax) == 0);
            
            fVecNew = obj.fVec(indfResMin:indfResMax);
            nPointsNew = length(fVecNew);
            bandwidthNew = fVecNew(end)-fVecNew(1);
            timeResolutionNew = 1/bandwidthNew;
            tVecNew = 0:timeResolutionNew:timeResolutionNew*(nPointsNew-1);
            % restrict CTF
            ctfNew = obj.CTF(:,indfResMin:indfResMax,:);
            
            % store to obj
            obj.fVec = fVecNew;
            obj.fMin = fResMin;
            obj.fMax = fResMax;
            obj.nPoints = nPointsNew;
            obj.CTF = ctfNew;
            obj.bandwidth = bandwidthNew;
            obj.timeResolution = timeResolutionNew;
            obj.tVec = tVecNew;
            
            obj.CIR = ctf2cir(obj);
        end
        
        function cir = ctf2cir(obj)
            if isempty(obj.CTF) == 0
                cir_tmp = zeros(size(obj.CTF));
                for i=1:obj.nTxAnt
                    for j=1:obj.nRxAnt
                        ctf_tmp = obj.CTF(i,:,j);
                        cir_tmp(i,:,j) = ifft(ctf_tmp);
                    end
                end
            else
                error('Only CIR is defined. Try the cir2ctf method.');
            end
            cir = cir_tmp;
        end
        
        function ctf = cir2ctf(obj)
            if isempty(obj.CIR) == 0
                ctf_tmp = zeros(size(obj.CIR));
                for i=1:obj.nTxAnt
                    for j=1:obj.nRxAnt
                        cir_tmp = obj.CIR(i,:,j);
                        ctf_tmp(i,:,j) = fft(cir_tmp);
                    end
                end
            else
                error('Only CTF is defined. Try the ctf2cir method.');
            end
            ctf = ctf_tmp;
        end
        
        function yLimits = getYlimits(obj)
            yLimits = zeros(1,2,2); % for CTF and CIR / finnaly in dB
            
            szCTF = size(obj.CTF);
            szCIR = size(obj.CIR);
            % get min/max value from CTF
            minCTFtmp = nan(szCTF(1), szCTF(3));
            maxCTFtmp = nan(szCTF(1), szCTF(3));
            for i=1:szCTF(1)
                for j = 1:szCTF(3)
                    minCTFtmp(i,j) = min(obj.CTF(i,:,j));
                    maxCTFtmp(i,j) = max(obj.CTF(i,:,j));
                end
            end
            % get min/max value from CIR
            minCIRtmp = nan(szCIR(1), szCIR(3));
            maxCIRtmp = nan(szCIR(1), szCIR(3));
            for i=1:szCIR(1)
                for j = 1:szCIR(3)
                    minCIRtmp(i,j) = min(obj.CIR(i,:,j));
                    maxCIRtmp(i,j) = max(obj.CIR(i,:,j));
                end
            end
            % save limits for CTF
            yLimitsTmp1 = [min(min(minCTFtmp)), max(max(maxCTFtmp))];
            yLimitsdBTmp1 = 20*log10(abs(yLimitsTmp1));
            yLimits(:,:,1) = [floor(yLimitsdBTmp1(1)/10)*10, ceil(yLimitsdBTmp1(2)/10)*10];
            % save limits for CIR
            yLimitsTmp2 = [min(min(minCIRtmp)), max(max(maxCIRtmp))];
            yLimitsdBTmp2 = 20*log10(abs(yLimitsTmp2));
            yLimits(:,:,2) = [floor(yLimitsdBTmp2(1)/10)*10, ceil(yLimitsdBTmp2(2)/10)*10];
        end
        
        function ctfResampled = resampleInterpCtf(obj, N)
            ctf = obj.CTF;
            
            nPoints = obj.nPoints;
            
            fMin = obj.fMin;
            fMax = obj.fMax;
            
            fVecNew = linspace(fMin, fMax, N);
            
            ctfNew = zeros(obj.nTxAnt,N,obj.nRxAnt);
            
            for i=1:obj.nTxAnt
                for j=1:obj.nRxAnt
                    ctfTemp = ctf(i,:,j);
%                     ctfNewTemp1 = interp1(fVec,ctfTemp,fVecNew,'linear');
%                     ctfNewTemp2 = interp1(fVec,ctfTemp,fVecNew,'pchip');
%                     ctfNewTemp3 = interp1(fVec,ctfTemp,fVecNew,'makima');
                    ctfNewTemp = interp1(obj.fVec,ctfTemp,fVecNew,'spline');
%                     plot()...
%                     legend('original','linear','pchip/cubic','makima','spline')
                    ctfNew(i,:,j) = ctfNewTemp;
                end
            end
            
            
        end
        
%         function [chanRmsPowMatrixF, chanRmsPowMatrixT] = getChanRmsPow(obj)
%            chanRmsPowMatrixF = zeros(obj.nTxAnt,obj.nRxAnt);
%            chanRmsPowMatrixT = zeros(obj.nTxAnt,obj.nRxAnt);
%            
%            for i=1:obj.nTxAnt
%                for j=1:obj.nRxAnt
%                     ctfTmp = obj.CTF(i,:,j);
%                     cirTmp = obj.CIR(i,:,j);
%                     
% %                     ctfRmsPowTemp = rms(ctfTmp);
%                     ctfRmsPowTemp = sqrt(sum(ctfTmp.*conj(ctfTmp))/length(ctfTmp));
% %                     cirRmsPowTemp = rms(cirTmp);
%                     cirRmsPowTemp = sqrt(sum(cirTmp.*conj(cirTmp))/length(cirTmp));
%                     
%                     chanRmsPowMatrixF(i,j) = ctfRmsPowTemp;
%                     chanRmsPowMatrixT(i,j) = cirRmsPowTemp;
%                end
%            end
%            
%         end
            
    end
end

% % Normalized (average) received power (NRP) (per sample) [-]
% P0 = (1/Nk)*sum(abs(magS).^2);
% 
% % Root mean square value
% P0_rms = sqrt(P0);
% % Normalized (average) received power (NRP) (per sample) [dB]
% P0_rms_dB = (1/2)*20*log10(P0); % in dB (1/2) -> sqrt;
