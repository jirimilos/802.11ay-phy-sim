function [coding] = RS_Blck_SPC_params(wifi_params, m_stbc, mcs_ind)
% Function for computing auxiliary Reed Solomon and Block coding parameters (LPSC PHY)
%
% Author:	Jiri Milos, DREL FEEC BUT, 2019
%

MCS = wifi_params.MCS(mcs_ind);
LENGTH = wifi_params.general.LENGTH;
N_cbps = MCS.N_cbps;

switch MCS.MCS
    case 'MCS25'
        coding.type{1} = 'RS';
        coding.type{2} = 'Blck';
        outer_code_params = [224 208];
        inner_code_params = [16 8];
    case 'MCS26'
        coding.type{1} = 'RS';
        coding.type{2} = 'Blck';
        outer_code_params = [224 208];
        inner_code_params = [12 8];
    case 'MCS27'
        coding.type{1} = 'RS';
        coding.type{2} = 'SPC';
        outer_code_params = [224 208];
        inner_code_params = [9 8];
    case 'MCS28'
        coding.type{1} = 'RS';
        coding.type{2} = 'Blck';
        outer_code_params = [224 208];
        inner_code_params = [16 8];
    case 'MCS29'
        coding.type{1} = 'RS';
        coding.type{2} = 'Blck';
        outer_code_params = [224 208];
        inner_code_params = [12 8];
    case 'MCS30'
        coding.type{1} = 'RS';
        coding.type{2} = 'SPC';
        outer_code_params = [224 208];
        inner_code_params = [9 8];
    case 'MCS31'
        coding.type{1} = 'RS';
        coding.type{2} = 'Blck';
        outer_code_params = [224 208];
        inner_code_params = [8 8];
    otherwise
        error('Wrong LPSC PHY MCS.');
end


% Generator matrices for inner block coding
switch inner_code_params(1)
    case 8
        G = eye(8);
    case 9
        G = [eye(8), ones(8, 1)];
    case 12
        rhs_part = [...
            1, 1, 0, 0;
            1, 0, 1, 0;
            0, 1, 0, 1;
            0, 0, 1, 1;
            1, 0, 1, 1;
            1, 1, 0, 1;
            1, 1, 1, 0;
            0, 1, 1, 1];
        G = [eye(8), rhs_part];
        clear rhs_part;
    case 16
        rhs_part = [...
            0, 1, 1, 0, 1, 0, 1, 0;
            0, 0, 1, 1, 0, 1, 0, 1;
            1, 0, 0, 1, 1, 0, 1, 0;
            0, 1, 0, 0, 1, 1, 0, 1;
            1, 0, 1, 0, 0, 1, 1, 0;
            0, 1, 0, 1, 0, 0, 1, 1;
            1, 0, 1, 0, 1, 0, 0, 1;
            1, 1, 0, 1, 0, 1, 0, 0];
        G = [eye(8), rhs_part];
        clear rhs_part;
    otherwise
        error('Wrong inner code parameter.');
end


% data is broken into blocks of length 208*8 bits

% Reed-Solomon outer code
% primitive polynomial
px_bin = [1 0 0 0 1 1 1 0 1]; % left msb
px_dec = bi2de(px_bin,'left-msb');
% number of bits per symbol
m = 8;
% maximum number of symbols in RS codeword within GF(2^m)
n_max = 2^m-1;
n1 = outer_code_params(1);
k1 = outer_code_params(2);
k1_max = n_max-(n1-k1);
% generator polynomial
gen_poly{1} = rsgenpoly(n_max, k1_max, px_dec, 0); % zero offset

% total number of Reed Solomon (RS) codewords
N_rs_cw = ceil(LENGTH/208);
% total number of RS encoded symbols
N_rs_encsym_tot = LENGTH+N_rs_cw*16;

rest_rs_block = LENGTH-(208*N_rs_cw);

num_rs_coders = 1;
if rest_rs_block < 0
    k2 = 208+rest_rs_block;
    n2 = k2+16;
    outer_code_params = [outer_code_params; n2, k2]; % RS coder parameters including the last codeword
    k2_max = n_max-(n2-k2);
    gen_poly{2} = rsgenpoly(n_max, k2_max, px_dec, 0); % zero offset
    num_rs_coders = 2;
end
% Block inner code
N = inner_code_params(1);
% total numbers of encoded bits
N_eb = N*(LENGTH+N_rs_cw*16);

% minimum number of blocks
N_blks_min = 18;
% TEST case (may be changed in future) - BRP is not used (see ch. 20.7.2.3.3.1)
BRP = 0;

% total number of 512 (in length) blocks
N_blks = ceil(N_eb/(N_cbps*392));
if (BRP ~= 0) && (N_blks < N_blks_min)
    N_blks = N_blks_min;
end

% Total number of zero block padding bits
N_blk_pad = (N_blks*N_cbps*392)-N_eb;

% Parity check matrix
[mg, ng] = size(G);
if mg ~= ng
    parmat = gen2par(G);
    % Syndrome decode table
    trt = syndtable(parmat);
else
    parmat = G;
    trt = [];
end

% save to structure
% outer code (Reed-Solomon)
coding.RS.N_cw = N_rs_cw;
coding.RS.N_encsym_tot = N_rs_encsym_tot;
coding.RS.n_k = outer_code_params;
coding.RS.n_max = n_max;
coding.RS.bitsPerSymbol = m;
coding.RS.P_bin = px_bin;
coding.RS.P_dec = px_dec;
coding.RS.G = gen_poly;
coding.RS.numEncoders = num_rs_coders;
% inner code (Block)
coding.Blck.N = N;
coding.Blck.G = G;
coding.Blck.ParCheckM = parmat;
coding.Blck.SyndTable = trt;
coding.Blck.N_eb = N_eb;
coding.Blck.n_k = inner_code_params;
coding.Blck.N_blks_min = N_blks_min;
coding.Blck.N_blks = N_blks;
coding.Blck.N_blk_pad = N_blk_pad;

end