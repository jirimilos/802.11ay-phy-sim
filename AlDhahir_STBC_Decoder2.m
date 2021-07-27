function [Y_out] = AlDhahir_STBC_Decoder2( SyPF, y_frames, h_mat, mode)

% "y_frames" consists of NTX received cplx frames without CP (NTX x SyPF)
% "h_mat" is a cplx channel impulse response matrix (NTX x SyPF)
% response in each ROW. (NTX x SyPF)

if(mode == 1)  % 2 TX and 1 RX antenna (2x1 MISO) Alamouti (Al-Dhahir)
    %   | s0  s1  |
    %   |-s1* s0* |
    % (Alamouti scheme)
    
    H0 = fft(h_mat(1,:));
    H1 = fft(h_mat(2,:));
    D0 = diag(H0 ,0);% create diagonal matrix
    D1 = diag(H1 ,0);
    Yx = zeros(SyPF, 2); % columns == FRAMES
    %     DD = abs(D1).^2 + abs(D2).^2; % diagonal matrix
    %     IDD = inv(DD);
    Y = fft(y_frames.').';
    % Aldhahir's Linear Processing on the RX side:
    %     for i=1:2:NF-1 % picking two Y-FRAMES at once
    aux0 = Y(1,:).';        % column vector
    aux1 = conj(Y(2,:)).';% column vector
    % RX Linear processing (STBC decoding):
    %     Yx(:,1) = IDD*(D1'*aux1 + D2*aux2);
    %     Yx(:,2) = IDD*(D2'*aux1 - D1*aux2);
    Yx(:,1) = (D0'*aux0 + D1*aux1);
    Yx(:,2) = (D1'*aux0 - D0*aux1);
    Y_out = Yx.'; % (2 x SyPF)
    
    
    
elseif(mode == 11)%  (2x2 MIMO) Alamouti (Al-Dhahir)
    %   | s0  s1  |
    %   |-s1* s0* |
    % (Alamouti scheme)
    H00 = fft(h_mat(1,1:SyPF));
    H01 = fft(h_mat(2,1:SyPF));
    H10 = fft(h_mat(3,1:SyPF));
    H11 = fft(h_mat(4,1:SyPF));
    
    D00 = diag(H00 ,0);% create diagonal matrix
    D01 = diag(H01 ,0);
    D10 = diag(H10 ,0);% create diagonal matrix
    D11 = diag(H11 ,0);
    
    Yx = zeros(SyPF, 2); % columns == FRAMES
    
    Y = fft(y_frames.').'; % (4 x SyPD)
    % Aldhahir's Linear Processing on the RX side:
    aux0_A = Y(1,:).';        % (SyPD x 1)
    aux1_A = conj(Y(2,:)).';% (SyPD x 1)
    aux0_B = Y(3,:).';        % (SyPD x 1)
    aux1_B = conj(Y(4,:)).';% (SyPD x 1)
    
    % RX Linear processing (STBC decoding):
    Yx(:,1) = (D00'*aux0_A + D10*aux1_A) + (D01'*aux0_B + D11*aux1_B);
    Yx(:,2) = (D10'*aux0_A - D00*aux1_A) + (D11'*aux0_B - D01*aux1_B);
    Y_out = Yx.'; % (2 x SyPF)
    
    
end