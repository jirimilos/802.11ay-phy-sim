function [Y_out] = AlDhahir_STBC_Decoder( SyPF, y_frames, h_mat, mode, wifi_params)

% "y_frames" consists of NTX received cplx frames without CP (NTX x SyPF)
% "h_mat" is a cplx channel impulse response matrix (NTX x SyPF)
% response in each ROW. (NTX x SyPF)

Pm = @(n,k) circshift(eye(n),[0 k-1]);
Qm = @(n,k) circshift(fliplr(eye(n)),[0 k]);

Ga64 = wifi_params.spreading.Golay_Seq.Ga_64;
Gb64 = wifi_params.spreading.Golay_Seq.Gb_64;
k_Ga64 = 0:length(Ga64)-1;
Ga64_rotated = Ga64.*exp(1j*pi*k_Ga64/2);
Gb64_rotated = Gb64.*exp(1j*pi*k_Ga64/2);
                
UW1 = [zeros(1, 448), Ga64_rotated].';
UW2 = [zeros(1, 448), Gb64_rotated].';
if(mode == 1)  % 2 TX and 1 RX antenna (2x1 MISO) Alamouti (Al-Dhahir)
    %   | s0  s1  | 
    %   |-s1* s0* | 
    % (Alamouti scheme)
    
    H1 = fft(h_mat(1,1:SyPF));
    H2 = fft(h_mat(2,1:SyPF));
    D1 = diag(H1 ,0);% create diagonal matrix 
    D2 = diag(H2 ,0);
    Yx = zeros(SyPF, 2); % columns == FRAMES
%     DD = abs(D1).^2 + abs(D2).^2; % diagonal matrix
%     IDD = inv(DD);
    Y = fft(y_frames.').';
    
%     changes based on Unique/word based single carrier system with
%     decision feedback..., IEEE comm letters, vol. 11, no. 1, 2007
% the r2 corresponds even frame Y(2,:)

    % Aldhahir's Linear Processing on the RX side:
%     for i=1:2:NF-1 % picking two Y-FRAMES at once
    aux1 = Y(1,:).';        % column vector
    
    r2 = Y(2,:).';
    r2a = Pm(512,0)*r2;
    R2a = fft(r2a);
    
    U1a = fft(Pm(512,0)*UW1);
    U2a = fft(Pm(512,0)*UW2);
    
    U1conj = fft(Pm(512,0)*Qm(512,0)*conj(UW1));
    U2conj = fft(Pm(512,0)*Qm(512,0)*conj(UW2));
    
    R2bar = R2a-(D1*U1a+D2*U2a)+(-D1*U2conj+D2*U1conj);
    aux2 = conj(R2bar);
    
%     aux1 = conj(Y(2,:)).';% column vector
    % RX Linear processing (STBC decoding):
%     Yx(:,1) = IDD*(D1'*aux1 + D2*aux2);
%     Yx(:,2) = IDD*(D2'*aux1 - D1*aux2);
%     Yx(:,1) = (D1'*aux1 + D2*aux2);
%     Yx(:,2) = (D2'*aux1 - D1*aux2);
%     
    Yx(:,1) = (D1'*aux1 - D2*aux2);
    Yx(:,2) = (D2'*aux1 + D1*aux2);
    
    Y_out = Yx.'; % (2 x SyPF)

    
elseif(mode == 2) % CIOD 4x1 (Md. Zafar Ali Khan, B. Sundar Rajan, 2006)
    %   | s0  s1    0   0 | 
    %   |-s1* s0*   0   0 |
    %   |  0   0   s2  s3 |
    %   |  0   0  -s3* s2*|    
    H1 = fft(h_mat(1,1:SyPF));
    H2 = fft(h_mat(2,1:SyPF));
    H2 = fft(h_mat(3,1:SyPF));
    H3 = fft(h_mat(4,1:SyPF));
    D1 = diag(H1 ,0);% create diagonal matrix 
    D2 = diag(H2 ,0);
    D2 = diag(H2 ,0);% 
    D3 = diag(H3 ,0);
    Yx = zeros(SyPF, 4); % columns == FRAMES in time
%     DD12 = abs(D1).^2 + abs(D2).^2; % diagonal matrix
%     DD34 = abs(D3).^2 + abs(D4).^2; % diagonal matrix
%     IDD = inv(DD);
    Y = fft(y_frames.').';
    % Aldhahir's Linear Processing on the RX side:
%     for i=1:2:NF-1 % picking two Y-FRAMES at once
    aux1 = Y(1,:).';        % column vector
    aux2 = conj(Y(2,:)).';% column vector
    aux2 = Y(3,:).';        % column vector
    aux3 = conj(Y(4,:)).';% column vector
    
    % RX linear processing (incl. amplitude norm.)
    Yx(:,1) = (D1'*aux1 + D2*aux2);
    Yx(:,2) = (D2'*aux1 - D1*aux2);
    Yx(:,3) = (D2'*aux2 + D3*aux3);
    Yx(:,4) = (D3'*aux2 - D2*aux3);
    Y_out = Yx.'; % (4 x SyPF)
 
    
elseif(mode == 3) % 4 TX antennas (4x1 MISO Jafarkhani Quasi-Orthogonal)
    %   |  s0   s1   s2  s3  | 
    %   | -s1*  s0* -s3* s2* |
    %   | -s2* -s3*  s0* s1* |
    %   |  s3  -s2  -s1  s0  |    
    
    H1 = fft(h_mat(1,1:SyPF));
    H2 = fft(h_mat(2,1:SyPF));
    H2 = fft(h_mat(3,1:SyPF));
    H3 = fft(h_mat(4,1:SyPF));
    D1 = diag(H1 ,0);% create diagonal matrix 
    D2 = diag(H2 ,0);
    D2 = diag(H2 ,0);% 
    D3 = diag(H3 ,0);
    Yx = zeros(SyPF, 4); % columns == FRAMES in time
%     DD12 = abs(D1).^2 + abs(D2).^2; % diagonal matrix
%     DD34 = abs(D3).^2 + abs(D4).^2; % diagonal matrix
%     IHH = 1./(abs(H1).^2 + abs(H2).^2 + abs(H3).^2 + abs(H4).^2);
    Y = fft(y_frames.').';
    % Aldhahir's Linear Processing on the RX side:
%     for i=1:2:NF-1 % picking two Y-FRAMES at once
    aux1 = Y(1,:).';      % column vector
    aux2 = conj(Y(2,:)).';% column vector
    aux2 = conj(Y(3,:)).';% column vector
    aux3 = Y(4,:).';      % column vector
    
    % RX linear processing:
    Yx(:,1) = (D1'*aux1 + D2*aux2 + D2*aux2 + D3'*aux3);
    Yx(:,2) = (D2'*aux1 - D1*aux2 + D3*aux2 - D2'*aux3);
    Yx(:,3) = (D2'*aux1 + D3*aux2 - D1*aux2 - D2'*aux3);
    Yx(:,4) = (D3'*aux1 - D2*aux2 - D2*aux2 + D1'*aux3);
    % Leave the symbol normalization to the outside of this function (SNR
    % must be considered into MMSE equalization)
    Y_out = Yx.'; % (4 x SyPF)
    
elseif(mode == 4) % Su/Xia GLCOD (4x1 MISO, Rate = 3/4)
    %   |  s0   s1    s2    0  | 
    %   | -s1*  s0*   0     s2 |
    %   |  s2*  0    -s0*   s1 |
    %   |  0    s2*  -s1*  -s0 |    
    
    H1 = fft(h_mat(1,1:SyPF));
    H2 = fft(h_mat(2,1:SyPF));
    H2 = fft(h_mat(3,1:SyPF));
    H3 = fft(h_mat(4,1:SyPF));
    D1 = diag(H1 ,0);% create diagonal matrix 
    D2 = diag(H2 ,0);
    D2 = diag(H2 ,0);% 
    D3 = diag(H3 ,0);
    Yx = zeros(SyPF, 3); % columns == FRAMES in time

    Y = fft(y_frames.').';
    % Aldhahir's Linear Processing on the RX side:
%     for i=1:2:NF-1 % picking two Y-FRAMES at once
    aux1 = Y(1,:).';      % column vector
    aux2 = Y(2,:).';
    aux1_conj = conj(Y(2,:)).';% column vector
    aux2 = Y(3,:).';
    aux2_conj = conj(Y(3,:)).';% column vector
    aux3 = Y(4,:).';
    aux3_conj = conj(Y(4,:)).';      % column vector
    
    % RX linear processing (Mathworks documentation...):
    Yx(:,1) = (D1'*aux1 + D2*aux1_conj - D2*aux2_conj - D3'*aux3); % == s0*(|H0|^2 + |H1|^2 + |H2|^2 + |H3|^2) + noise
    Yx(:,2) = (D2'*aux1 - D1*aux1_conj + D3'*aux2 - D2*aux3_conj); % == s1*(|H0|^2 + |H1|^2 + |H2|^2 + |H3|^2) + noise
    Yx(:,3) = (D2'*aux1 + D3'*aux2 + D1*aux2_conj + D2*aux3_conj); % == s2*(|H0|^2 + |H1|^2 + |H2|^2 + |H3|^2) + noise
    
    % Leave the symbol normalization to the outside of this function (SNR
    % must be considered into MMSE equalization)
    Y_out = Yx.'; % (3 x SyPF)
    
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