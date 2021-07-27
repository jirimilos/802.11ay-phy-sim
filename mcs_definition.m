%% Script defines Modulation and Coding Schemes
%
% Author:	Jiri Milos, DREL FEEC BUT, 2018
%
global wifi_params

switch wifi_standard
%% IEEE 802.11ad     
    case '802dot11ad'
%             case 'C-PHY'
                wifi_params.MCS(1).phy_type = 'Ctrl'; % stands for Control
                wifi_params.MCS(1).N_ss = 1;
                wifi_params.MCS(1).field_base = 0; % Base MCS field
                wifi_params.MCS(1).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(1).MCS = 'MCS0';
                wifi_params.MCS(1).CR = 1/2;
                wifi_params.MCS(1).Repetition = 1;
                wifi_params.MCS(1).M = 2; % pi/2-DBPSK
                wifi_params.MCS(1).coding_type = 'LDPC';
                wifi_params.MCS(1).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(1).N_cbpb = 2*168; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(1).DataRate = 27.5e6; % Max. theoretical throughput, or data rate
%             case 'SC-PHY'
                % Single Carrier
                % MCS 1 -------------------------------------------------------------------
                wifi_params.MCS(2).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(2).N_ss = 1;
                wifi_params.MCS(2).field_base = 1; % Base MCS field
                wifi_params.MCS(2).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(2).MCS = 'MCS1';
                wifi_params.MCS(2).CR = 1/2;
                wifi_params.MCS(2).Repetition = 2;
                wifi_params.MCS(2).M = 2; % pi/2-BPSK
                wifi_params.MCS(2).coding_type = 'LDPC';
                wifi_params.MCS(2).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(2).N_cbpb = 448; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(2).DataRate = 385e6; % Max. theoretical throughput, or data rate
                % MCS 2 -------------------------------------------------------------------
                wifi_params.MCS(3).phy_type = 'SC';
                wifi_params.MCS(3).N_ss = 1;
                wifi_params.MCS(3).field_base = 2; % Base MCS field
                wifi_params.MCS(3).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(3).MCS = 'MCS2';
                wifi_params.MCS(3).CR = 1/2;
                wifi_params.MCS(3).Repetition = 1;
                wifi_params.MCS(3).M = 2; % pi/2-BPSK
                wifi_params.MCS(3).coding_type = 'LDPC';
                wifi_params.MCS(3).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(3).N_cbpb = 448; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(3).DataRate = 770e6; % Max. theoretical throughput, or data rate
                % MCS 3 -------------------------------------------------------------------
                wifi_params.MCS(4).phy_type = 'SC';
                wifi_params.MCS(4).N_ss = 1;
                wifi_params.MCS(4).field_base = 3; % Base MCS field
                wifi_params.MCS(4).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(4).MCS = 'MCS3';
                wifi_params.MCS(4).CR = 5/8;
                wifi_params.MCS(4).Repetition = 1;
                wifi_params.MCS(4).M = 2; % pi/2-BPSK
                wifi_params.MCS(4).coding_type = 'LDPC';
                wifi_params.MCS(4).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(4).N_cbpb = 448; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(4).DataRate = 962.5e6; % Max. theoretical throughput, or data rate
                % MCS 4 -------------------------------------------------------------------
                wifi_params.MCS(5).phy_type = 'SC';
                wifi_params.MCS(5).N_ss = 1;
                wifi_params.MCS(5).field_base = 4; % Base MCS field
                wifi_params.MCS(5).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(5).MCS = 'MCS4';
                wifi_params.MCS(5).CR = 3/4;
                wifi_params.MCS(5).Repetition = 1;
                wifi_params.MCS(5).M = 2; % pi/2-BPSK
                wifi_params.MCS(5).coding_type = 'LDPC';
                wifi_params.MCS(5).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(5).N_cbpb = 448; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(5).DataRate = 1155e6; % Max. theoretical throughput, or data rate
                % MCS 5 -------------------------------------------------------------------
                wifi_params.MCS(6).phy_type = 'SC';
                wifi_params.MCS(6).N_ss = 1;
                wifi_params.MCS(6).field_base = 5; % Base MCS field
                wifi_params.MCS(6).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(6).MCS = 'MCS5';
                wifi_params.MCS(6).CR = 13/16;
                wifi_params.MCS(6).Repetition = 1;
                wifi_params.MCS(6).M = 2; % pi/2-BPSK
                wifi_params.MCS(6).coding_type = 'LDPC';
                wifi_params.MCS(6).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(6).N_cbpb = 448; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(6).DataRate = 1251.25e6; % Max. theoretical throughput, or data rate
                % MCS 6 -------------------------------------------------------------------
                wifi_params.MCS(7).phy_type = 'SC';
                wifi_params.MCS(7).N_ss = 1;
                wifi_params.MCS(7).field_base = 6; % Base MCS field
                wifi_params.MCS(7).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(7).MCS = 'MCS6';
                wifi_params.MCS(7).CR = 1/2;
                wifi_params.MCS(7).Repetition = 1;
                wifi_params.MCS(7).M = 4; % pi/2-QPSK
                wifi_params.MCS(7).coding_type = 'LDPC';
                wifi_params.MCS(7).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(7).N_cbpb = 896; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(7).DataRate = 1540e6; % Max. theoretical throughput, or data rate
                % MCS 7 -------------------------------------------------------------------
                wifi_params.MCS(8).phy_type = 'SC';
                wifi_params.MCS(8).N_ss = 1;
                wifi_params.MCS(8).field_base = 7; % Base MCS field
                wifi_params.MCS(8).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(8).MCS = 'MCS7';
                wifi_params.MCS(8).CR = 5/8;
                wifi_params.MCS(8).Repetition = 1;
                wifi_params.MCS(8).M = 4; % pi/2-QPSK
                wifi_params.MCS(8).coding_type = 'LDPC';
                wifi_params.MCS(8).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(8).N_cbpb = 896; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(8).DataRate = 1925e6; % Max. theoretical throughput, or data rate
                % MCS 8 -------------------------------------------------------------------
                wifi_params.MCS(9).phy_type = 'SC';
                wifi_params.MCS(9).N_ss = 1;
                wifi_params.MCS(9).field_base = 8; % Base MCS field
                wifi_params.MCS(9).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(9).MCS = 'MCS8';
                wifi_params.MCS(9).CR = 3/4;
                wifi_params.MCS(9).Repetition = 1;
                wifi_params.MCS(9).M = 4; % pi/2-QPSK
                wifi_params.MCS(9).coding_type = 'LDPC';
                wifi_params.MCS(9).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(9).N_cbpb = 896; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(9).DataRate = 2310e6; % Max. theoretical throughput, or data rate
                % MCS 9 -------------------------------------------------------------------
                wifi_params.MCS(10).phy_type = 'SC';
                wifi_params.MCS(10).N_ss = 1;
                wifi_params.MCS(10).field_base = 9; % Base MCS field
                wifi_params.MCS(10).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(10).MCS = 'MCS9';
                wifi_params.MCS(10).CR = 13/16;
                wifi_params.MCS(10).Repetition = 1;
                wifi_params.MCS(10).M = 4; % pi/2-QPSK
                wifi_params.MCS(10).coding_type = 'LDPC';
                wifi_params.MCS(10).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(10).N_cbpb = 896; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(10).DataRate = 2502.5e6; % Max. theoretical throughput, or data rate
                % MCS 9.1 -----------------------------------------------------------------
                wifi_params.MCS(11).phy_type = 'SC';
                wifi_params.MCS(11).N_ss = 1;
                wifi_params.MCS(11).field_base = 6; % Base MCS field
                wifi_params.MCS(11).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(11).MCS = 'MCS9.1';
                wifi_params.MCS(11).CR = 7/8;
                wifi_params.MCS(11).Repetition = 1;
                wifi_params.MCS(11).M = 4; % pi/2-QPSK
                wifi_params.MCS(11).coding_type = 'LDPC';
                wifi_params.MCS(11).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(11).N_cbpb = 896; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(11).DataRate = 2695e6; % Max. theoretical throughput, or data rate
                % MCS 10 ------------------------------------------------------------------
                wifi_params.MCS(12).phy_type = 'SC';
                wifi_params.MCS(12).N_ss = 1;
                wifi_params.MCS(12).field_base = 10; % Base MCS field
                wifi_params.MCS(12).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(12).MCS = 'MCS10';
                wifi_params.MCS(12).CR = 1/2;
                wifi_params.MCS(12).Repetition = 1;
                wifi_params.MCS(12).M = 16; % pi/2-16QAM
                wifi_params.MCS(12).coding_type = 'LDPC';
                wifi_params.MCS(12).N_cbps = 4; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(12).N_cbpb = 1792; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(12).DataRate = 3080e6; % Max. theoretical throughput, or data rate
                % MCS 11 ------------------------------------------------------------------
                wifi_params.MCS(13).phy_type = 'SC';
                wifi_params.MCS(13).N_ss = 1;
                wifi_params.MCS(13).field_base = 11; % Base MCS field
                wifi_params.MCS(13).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(13).MCS = 'MCS11';
                wifi_params.MCS(13).CR = 5/8;
                wifi_params.MCS(13).Repetition = 1;
                wifi_params.MCS(13).M = 16; % pi/2-16QAM
                wifi_params.MCS(13).coding_type = 'LDPC';
                wifi_params.MCS(13).N_cbps = 4; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(13).N_cbpb = 1792; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(13).DataRate = 3850e6; % Max. theoretical throughput, or data rate
                % MCS 12 ------------------------------------------------------------------
                wifi_params.MCS(14).phy_type = 'SC';
                wifi_params.MCS(14).N_ss = 1;
                wifi_params.MCS(14).field_base = 12; % Base MCS field
                wifi_params.MCS(14).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(14).MCS = 'MCS12';
                wifi_params.MCS(14).CR = 3/4;
                wifi_params.MCS(14).Repetition = 1;
                wifi_params.MCS(14).M = 16; % pi/2-16QAM
                wifi_params.MCS(14).coding_type = 'LDPC';
                wifi_params.MCS(14).N_cbps = 4; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(14).N_cbpb = 1792; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(14).DataRate = 4620e6; % Max. theoretical throughput, or data rate
                % MCS 12.1 ----------------------------------------------------------------
                wifi_params.MCS(15).phy_type = 'SC';
                wifi_params.MCS(15).N_ss = 1;
                wifi_params.MCS(15).field_base = 7; % Base MCS field
                wifi_params.MCS(15).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(15).MCS = 'MCS12.1';
                wifi_params.MCS(15).CR = 13/16;
                wifi_params.MCS(15).Repetition = 1;
                wifi_params.MCS(15).M = 16; % pi/2-16QAM
                wifi_params.MCS(15).coding_type = 'LDPC';
                wifi_params.MCS(15).N_cbps = 4; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(15).N_cbpb = 1792; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(15).DataRate = 5005e6; % Max. theoretical throughput, or data rate
                % MCS 12.2 ----------------------------------------------------------------
                wifi_params.MCS(16).phy_type = 'SC';
                wifi_params.MCS(16).N_ss = 1;
                wifi_params.MCS(16).field_base = 8; % Base MCS field
                wifi_params.MCS(16).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(16).MCS = 'MCS12.2';
                wifi_params.MCS(16).CR = 7/8;
                wifi_params.MCS(16).Repetition = 1;
                wifi_params.MCS(16).M = 16; % pi/2-16QAM
                wifi_params.MCS(16).coding_type = 'LDPC';
                wifi_params.MCS(16).N_cbps = 4; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(16).N_cbpb = 1792; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(16).DataRate = 5390e6; % Max. theoretical throughput, or data rate
                % MCS 12.3 ----------------------------------------------------------------
                wifi_params.MCS(17).phy_type = 'SC';
                wifi_params.MCS(17).N_ss = 1;
                wifi_params.MCS(17).field_base = 9; % Base MCS field
                wifi_params.MCS(17).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(17).MCS = 'MCS12.3';
                wifi_params.MCS(17).CR = 5/8;
                wifi_params.MCS(17).Repetition = 1;
                wifi_params.MCS(17).M = 64; % pi/2-64QAM
                wifi_params.MCS(17).coding_type = 'LDPC';
                wifi_params.MCS(17).N_cbps = 6; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(17).N_cbpb = 2688; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(17).DataRate = 5775e6; % Max. theoretical throughput, or data rate
                % MCS 12.4 ----------------------------------------------------------------
                wifi_params.MCS(18).phy_type = 'SC';
                wifi_params.MCS(18).N_ss = 1;
                wifi_params.MCS(18).field_base = 10; % Base MCS field
                wifi_params.MCS(18).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(18).MCS = 'MCS12.4';
                wifi_params.MCS(18).CR = 3/4;
                wifi_params.MCS(18).Repetition = 1;
                wifi_params.MCS(18).M = 64; % pi/2-64QAM
                wifi_params.MCS(18).coding_type = 'LDPC';
                wifi_params.MCS(18).N_cbps = 6; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(18).N_cbpb = 2688; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(18).DataRate = 6390e6; % Max. theoretical throughput, or data rate
                % MCS 12.5 ----------------------------------------------------------------
                wifi_params.MCS(19).phy_type = 'SC';
                wifi_params.MCS(19).N_ss = 1;
                wifi_params.MCS(19).field_base = 11; % Base MCS field
                wifi_params.MCS(19).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(19).MCS = 'MCS12.5';
                wifi_params.MCS(19).CR = 13/16;
                wifi_params.MCS(19).Repetition = 1;
                wifi_params.MCS(19).M = 64; % pi/2-64QAM
                wifi_params.MCS(19).coding_type = 'LDPC';
                wifi_params.MCS(19).N_cbps = 6; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(19).N_cbpb = 2688; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(19).DataRate = 7507.5e6; % Max. theoretical throughput, or data rate
                % MCS 12.6 ----------------------------------------------------------------
                wifi_params.MCS(20).phy_type = 'SC';
                wifi_params.MCS(20).N_ss = 1;
                wifi_params.MCS(20).field_base = 12; % Base MCS field
                wifi_params.MCS(20).field_extended = 1; % Extended SC MCS indication field
                wifi_params.MCS(20).MCS = 'MCS12.6';
                wifi_params.MCS(20).CR = 7/8;
                wifi_params.MCS(20).Repetition = 1;
                wifi_params.MCS(20).M = 64; % pi/2-64QAM
                wifi_params.MCS(20).coding_type = 'LDPC';
                wifi_params.MCS(20).N_cbps = 6; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(20).N_cbpb = 2688; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(20).DataRate = 8085e6; % Max. theoretical throughput, or data rate
                % case 'OFDM'
                % MCS 13 ------------------------------------------------------------------
                wifi_params.MCS(21).phy_type = 'OFDM';
                wifi_params.MCS(21).N_ss = 1;
                wifi_params.MCS(21).field_base = 13; % Base MCS field
                wifi_params.MCS(21).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(21).MCS = 'MCS13';
                wifi_params.MCS(21).CR = 1/2;
                wifi_params.MCS(21).Repetition = 1;
                wifi_params.MCS(21).M = 4; % SQPSK
                wifi_params.MCS(21).coding_type = 'LDPC';
                wifi_params.MCS(21).N_cbps = 336; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(21).N_dbps = 168; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(21).N_bpsc = 1; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(21).DataRate = 693e6; % Max. theoretical throughput, or data rate
                % MCS 14 ------------------------------------------------------------------
                wifi_params.MCS(22).phy_type = 'OFDM';
                wifi_params.MCS(22).N_ss = 1;
                wifi_params.MCS(22).field_base = 14; % Base MCS field
                wifi_params.MCS(22).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(22).MCS = 'MCS14';
                wifi_params.MCS(22).CR = 5/8;
                wifi_params.MCS(22).Repetition = 1;
                wifi_params.MCS(22).M = 4; % SQPSK
                wifi_params.MCS(22).coding_type = 'LDPC';
                wifi_params.MCS(22).N_cbps = 336; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(22).N_dbps = 210; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(22).N_bpsc = 1; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(22).DataRate = 866.25e6; % Max. theoretical throughput, or data rate
                % MCS 15 ------------------------------------------------------------------
                wifi_params.MCS(23).phy_type = 'OFDM';
                wifi_params.MCS(23).N_ss = 1;
                wifi_params.MCS(23).field_base = 15; % Base MCS field
                wifi_params.MCS(23).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(23).MCS = 'MCS15';
                wifi_params.MCS(23).CR = 1/2;
                wifi_params.MCS(23).Repetition = 1;
                wifi_params.MCS(23).M = 4; % QPSK
                wifi_params.MCS(23).coding_type = 'LDPC';
                wifi_params.MCS(23).N_cbps = 672; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(23).N_dbps = 336; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(23).N_bpsc = 2; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(23).DataRate = 1386e6; % Max. theoretical throughput, or data rate
                % MCS 16 ------------------------------------------------------------------
                wifi_params.MCS(24).phy_type = 'OFDM';
                wifi_params.MCS(24).N_ss = 1;
                wifi_params.MCS(24).field_base = 16; % Base MCS field
                wifi_params.MCS(24).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(24).MCS = 'MCS16';
                wifi_params.MCS(24).CR = 5/8;
                wifi_params.MCS(24).Repetition = 1;
                wifi_params.MCS(24).M = 4; % QPSK
                wifi_params.MCS(24).coding_type = 'LDPC';
                wifi_params.MCS(24).N_cbps = 672; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(24).N_dbps = 420; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(24).N_bpsc = 2; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(24).DataRate = 1732.5e6; % Max. theoretical throughput, or data rate
                % MCS 17 ------------------------------------------------------------------
                wifi_params.MCS(25).phy_type = 'OFDM';
                wifi_params.MCS(25).N_ss = 1;
                wifi_params.MCS(25).field_base = 17; % Base MCS field
                wifi_params.MCS(25).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(25).MCS = 'MCS17';
                wifi_params.MCS(25).CR = 3/4;
                wifi_params.MCS(25).Repetition = 1;
                wifi_params.MCS(25).M = 4; % QPSK
                wifi_params.MCS(25).coding_type = 'LDPC';
                wifi_params.MCS(25).N_cbps = 672; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(25).N_dbps = 504; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(25).N_bpsc = 2; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(25).DataRate = 2079e6; % Max. theoretical throughput, or data rate
                % MCS 18 ------------------------------------------------------------------
                wifi_params.MCS(26).phy_type = 'OFDM';
                wifi_params.MCS(26).N_ss = 1;
                wifi_params.MCS(26).field_base = 18; % Base MCS field
                wifi_params.MCS(26).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(26).MCS = 'MCS18';
                wifi_params.MCS(26).CR = 1/2;
                wifi_params.MCS(26).Repetition = 1;
                wifi_params.MCS(26).M = 16; % 16QAM
                wifi_params.MCS(26).coding_type = 'LDPC';
                wifi_params.MCS(26).N_cbps = 1344; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(26).N_dbps = 672; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(26).N_bpsc = 4; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(26).DataRate = 2772e6; % Max. theoretical throughput, or data rate
                % MCS 19 ------------------------------------------------------------------
                wifi_params.MCS(27).phy_type = 'OFDM';
                wifi_params.MCS(27).N_ss = 1;
                wifi_params.MCS(27).field_base = 19; % Base MCS field
                wifi_params.MCS(27).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(27).MCS = 'MCS19';
                wifi_params.MCS(27).CR = 5/8;
                wifi_params.MCS(27).Repetition = 1;
                wifi_params.MCS(27).M = 16; % 16QAM
                wifi_params.MCS(27).coding_type = 'LDPC';
                wifi_params.MCS(27).N_cbps = 1344; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(27).N_dbps = 840; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(27).N_bpsc = 4; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(27).DataRate = 3465e6; % Max. theoretical throughput, or data rate
                % MCS 20 ------------------------------------------------------------------
                wifi_params.MCS(28).phy_type = 'OFDM';
                wifi_params.MCS(28).N_ss = 1;
                wifi_params.MCS(28).field_base = 20; % Base MCS field
                wifi_params.MCS(28).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(28).MCS = 'MCS20';
                wifi_params.MCS(28).CR = 3/4;
                wifi_params.MCS(28).Repetition = 1;
                wifi_params.MCS(28).M = 16; % 16QAM
                wifi_params.MCS(28).coding_type = 'LDPC';
                wifi_params.MCS(28).N_cbps = 1344; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(28).N_dbps = 1008; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(28).N_bpsc = 4; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(28).DataRate = 4158e6; % Max. theoretical throughput, or data rate
                % MCS 21 ------------------------------------------------------------------
                wifi_params.MCS(29).phy_type = 'OFDM';
                wifi_params.MCS(29).N_ss = 1;
                wifi_params.MCS(29).field_base = 21; % Base MCS field
                wifi_params.MCS(29).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(29).MCS = 'MCS21';
                wifi_params.MCS(29).CR = 13/16;
                wifi_params.MCS(29).Repetition = 1;
                wifi_params.MCS(29).M = 16; % 16QAM
                wifi_params.MCS(29).coding_type = 'LDPC';
                wifi_params.MCS(29).N_cbps = 1344; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(29).N_dbps = 1092; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(29).N_bpsc = 4; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(29).DataRate = 4504.5e6; % Max. theoretical throughput, or data rate
                % MCS 22 ------------------------------------------------------------------
                wifi_params.MCS(30).phy_type = 'OFDM';
                wifi_params.MCS(30).N_ss = 1;
                wifi_params.MCS(30).field_base = 22; % Base MCS field
                wifi_params.MCS(30).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(30).MCS = 'MCS22';
                wifi_params.MCS(30).CR = 5/8;
                wifi_params.MCS(30).Repetition = 1;
                wifi_params.MCS(30).M = 64; % 64QAM
                wifi_params.MCS(30).coding_type = 'LDPC';
                wifi_params.MCS(30).N_cbps = 2016; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(30).N_dbps = 1260; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(30).N_bpsc = 6; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(30).DataRate = 5197.5e6; % Max. theoretical throughput, or data rate
                % MCS 23 ------------------------------------------------------------------
                wifi_params.MCS(31).phy_type = 'OFDM';
                wifi_params.MCS(31).N_ss = 1;
                wifi_params.MCS(31).field_base = 23; % Base MCS field
                wifi_params.MCS(31).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(31).MCS = 'MCS23';
                wifi_params.MCS(31).CR = 3/4;
                wifi_params.MCS(31).Repetition = 1;
                wifi_params.MCS(31).M = 64; % 64QAM
                wifi_params.MCS(31).coding_type = 'LDPC';
                wifi_params.MCS(31).N_cbps = 2016; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(31).N_dbps = 1512; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(31).N_bpsc = 6; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(31).DataRate = 6237e6; % Max. theoretical throughput, or data rate
                % MCS 24 ------------------------------------------------------------------
                wifi_params.MCS(32).phy_type = 'OFDM';
                wifi_params.MCS(32).N_ss = 1;
                wifi_params.MCS(32).field_base = 24; % Base MCS field
                wifi_params.MCS(32).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(32).MCS = 'MCS24';
                wifi_params.MCS(32).CR = 13/16;
                wifi_params.MCS(32).Repetition = 1;
                wifi_params.MCS(32).M = 64; % 64QAM
                wifi_params.MCS(32).coding_type = 'LDPC';
                wifi_params.MCS(32).N_cbps = 2016; % Number of coded bits per symbol, see Table 20-14                
                wifi_params.MCS(32).N_dbps = 1638; % Number of data bits per symbol, see Table 20-14                
                wifi_params.MCS(32).N_bpsc = 6; % Number of bits per subcarrier, see Table 20-14
                wifi_params.MCS(32).DataRate = 6756.75e6; % Max. theoretical throughput, or data rate
                % case 'LP-SC-PHY'
                % MCS 25 ----------------------------------------------------------------
                wifi_params.MCS(33).phy_type = 'LPSC';
                wifi_params.MCS(33).N_ss = 1;
                wifi_params.MCS(33).field_base = 25; % Base MCS field
                wifi_params.MCS(33).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(33).MCS = 'MCS25';
                wifi_params.MCS(33).CR = 13/28;
                wifi_params.MCS(33).Repetition = 1;
                wifi_params.MCS(33).M = 2; % pi/2-BPSK
                wifi_params.MCS(33).coding_type = 'RS+Blck';
                wifi_params.MCS(33).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(33).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(33).DataRate = 626e6; % Max. theoretical throughput, or data rate
                % MCS 26 ----------------------------------------------------------------
                wifi_params.MCS(34).phy_type = 'LPSC';
                wifi_params.MCS(34).N_ss = 1;
                wifi_params.MCS(34).field_base = 26; % Base MCS field
                wifi_params.MCS(34).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(34).MCS = 'MCS26';
                wifi_params.MCS(34).CR = 13/21;
                wifi_params.MCS(34).Repetition = 1;
                wifi_params.MCS(34).M = 2; % pi/2-BPSK
                wifi_params.MCS(34).coding_type = 'RS+Blck';
                wifi_params.MCS(34).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(34).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(34).DataRate = 834e6; % Max. theoretical throughput, or data rate
                % MCS 27 ----------------------------------------------------------------
                wifi_params.MCS(35).phy_type = 'LPSC';
                wifi_params.MCS(35).N_ss = 1;
                wifi_params.MCS(35).field_base = 27; % Base MCS field
                wifi_params.MCS(35).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(35).MCS = 'MCS27';
                wifi_params.MCS(35).CR = 52/63;
                wifi_params.MCS(35).Repetition = 1;
                wifi_params.MCS(35).M = 2; % pi/2-BPSK
                wifi_params.MCS(35).coding_type = 'RS+SPC';
                wifi_params.MCS(35).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(35).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(35).DataRate = 1112e6; % Max. theoretical throughput, or data rate
                % MCS 28 ----------------------------------------------------------------
                wifi_params.MCS(36).phy_type = 'LPSC';
                wifi_params.MCS(36).N_ss = 1;
                wifi_params.MCS(36).field_base = 28; % Base MCS field
                wifi_params.MCS(36).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(36).MCS = 'MCS28';
                wifi_params.MCS(36).CR = 13/28;
                wifi_params.MCS(36).Repetition = 1;
                wifi_params.MCS(36).M = 4; % pi/2-QPSK
                wifi_params.MCS(36).coding_type = 'RS+Blck';
                wifi_params.MCS(36).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(36).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(36).DataRate = 1251e6; % Max. theoretical throughput, or data rate
                % MCS 29 ----------------------------------------------------------------
                wifi_params.MCS(37).phy_type = 'LPSC';
                wifi_params.MCS(37).N_ss = 1;
                wifi_params.MCS(37).field_base = 29; % Base MCS field
                wifi_params.MCS(37).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(37).MCS = 'MCS29';
                wifi_params.MCS(37).CR = 13/21;
                wifi_params.MCS(37).Repetition = 1;
                wifi_params.MCS(37).M = 4; % pi/2-QPSK
                wifi_params.MCS(37).coding_type = 'RS+Blck';
                wifi_params.MCS(37).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(37).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(37).DataRate = 1668e6; % Max. theoretical throughput, or data rate
                % MCS 30 ----------------------------------------------------------------
                wifi_params.MCS(38).phy_type = 'LPSC';
                wifi_params.MCS(38).N_ss = 1;
                wifi_params.MCS(38).field_base = 30; % Base MCS field
                wifi_params.MCS(38).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(38).MCS = 'MCS30';
                wifi_params.MCS(38).CR = 52/63;
                wifi_params.MCS(38).Repetition = 1;
                wifi_params.MCS(38).M = 4; % pi/2-QPSK
                wifi_params.MCS(38).coding_type = 'RS+SPC';
                wifi_params.MCS(38).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(38).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21                
                wifi_params.MCS(38).DataRate = 2224e6; % Max. theoretical throughput, or data rate
                % MCS 31 ----------------------------------------------------------------
                wifi_params.MCS(39).phy_type = 'LPSC';
                wifi_params.MCS(39).N_ss = 1;
                wifi_params.MCS(39).field_base = 31; % Base MCS field
                wifi_params.MCS(39).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(39).MCS = 'MCS31';
                wifi_params.MCS(39).CR = 13/14;
                wifi_params.MCS(39).Repetition = 1;
                wifi_params.MCS(39).M = 4; % pi/2-QPSK
                wifi_params.MCS(39).coding_type = 'RS+Blck';
                wifi_params.MCS(39).N_cbps = 2; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(39).N_cbpb = 392; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(39).DataRate = 2503e6; % Max. theoretical throughput, or data rate
%% IEEE 802.11ay
    case '802dot11ay'
%             case 'C-PHY'
                wifi_params.MCS(1).phy_type = 'Ctrl'; % stands for Control
                wifi_params.MCS(1).N_ss = 1; % Possible numbers of spatial stream
                wifi_params.MCS(1).N_CB = 1; % Possible numbers of bonded channels
                wifi_params.MCS(1).field_base = 0; % Base MCS field
                wifi_params.MCS(1).field_extended = 0; % Extended SC MCS indication field
                wifi_params.MCS(1).MCS = 'MCS0';
                wifi_params.MCS(1).CR = 1/2;
                wifi_params.MCS(1).Repetition = 1;
                wifi_params.MCS(1).M = 2; % pi/2-DBPSK
                wifi_params.MCS(1).coding_type = 'LDPC';
                wifi_params.MCS(1).N_cbps = 1; % Number of coded bits per symbol, see Table 20-19
                wifi_params.MCS(1).N_cbpb = 2*168; % Number of coded bits per symbol block, see Table 20-21
                wifi_params.MCS(1).DataRate = 27.5*[412.5/385, 1, 330/385]*1e6; % TBD, Max. theoretical throughput, or data rate (for short, normal or long GI length)
%             case 'SC-PHY'
                % Single Carrier
                % MCS 1 -------------------------------------------------------------------
                wifi_params.MCS(2).phy_type = 'SC'; % stands for Single carrier
                wifi_params.MCS(2).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(2).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(2).field_base = 1; % Base MCS field
                wifi_params.MCS(2).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(2).MCS = 'MCS1';
                wifi_params.MCS(2).CR = 1/2;
                wifi_params.MCS(2).Repetition = 2;
                wifi_params.MCS(2).M = 2; % pi/2-BPSK
                wifi_params.MCS(2).coding_type = 'LDPC';
                wifi_params.MCS(2).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(2).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(2).DataRate = [412.5, 385.0, 330.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 2 -------------------------------------------------------------------
                wifi_params.MCS(3).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(3).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(3).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(3).field_base = 2; % Base MCS field
                wifi_params.MCS(3).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(3).MCS = 'MCS2';
                wifi_params.MCS(3).CR = 1/2;
                wifi_params.MCS(3).Repetition = 1;
                wifi_params.MCS(3).M = 2; % pi/2-BPSK
                wifi_params.MCS(3).coding_type = 'LDPC';
                wifi_params.MCS(3).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(3).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(3).DataRate = [825.0, 770.0, 660.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 3 -------------------------------------------------------------------
                wifi_params.MCS(4).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(4).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(4).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(4).field_base = 3; % Base MCS field
                wifi_params.MCS(4).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(4).MCS = 'MCS3';
                wifi_params.MCS(4).CR = 5/8;
                wifi_params.MCS(4).Repetition = 1;
                wifi_params.MCS(4).M = 2; % pi/2-BPSK
                wifi_params.MCS(4).coding_type = 'LDPC';
                wifi_params.MCS(4).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(4).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(4).DataRate = [1031.25, 962.5, 825.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 4 -------------------------------------------------------------------
                wifi_params.MCS(5).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(5).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(5).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(5).field_base = 4; % Base MCS field
                wifi_params.MCS(5).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(5).MCS = 'MCS4';
                wifi_params.MCS(5).CR = 3/4;
                wifi_params.MCS(5).Repetition = 1;
                wifi_params.MCS(5).M = 2; % pi/2-BPSK
                wifi_params.MCS(5).coding_type = 'LDPC';
                wifi_params.MCS(5).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(5).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(5).DataRate = [1237.5, 1155.0, 990.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 5 -------------------------------------------------------------------
                wifi_params.MCS(6).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(6).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(6).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(6).field_base = 5; % Base MCS field
                wifi_params.MCS(6).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(6).MCS = 'MCS5';
                wifi_params.MCS(6).CR = 13/16;
                wifi_params.MCS(6).Repetition = 1;
                wifi_params.MCS(6).M = 2; % pi/2-BPSK
                wifi_params.MCS(6).coding_type = 'LDPC';
                wifi_params.MCS(6).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(6).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(6).DataRate = [1340.63, 1251.25, 1072.5]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 6 -------------------------------------------------------------------
                wifi_params.MCS(7).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(7).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(7).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(7).field_base = 6; % Base MCS field
                wifi_params.MCS(7).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(7).MCS = 'MCS6';
                wifi_params.MCS(7).CR = 7/8;
                wifi_params.MCS(7).Repetition = 1;
                wifi_params.MCS(7).M = 2; % pi/2-BPSK
                wifi_params.MCS(7).coding_type = 'LDPC';
                wifi_params.MCS(7).N_cbps = 1; % Number of coded bits per symbol
                wifi_params.MCS(7).N_cbpb = 1*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(7).DataRate = [1443.75, 1347.5, 1155.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 7 -------------------------------------------------------------------
                wifi_params.MCS(8).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(8).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(8).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(8).field_base = 7; % Base MCS field
                wifi_params.MCS(8).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(8).MCS = 'MCS7';
                wifi_params.MCS(8).CR = 1/2;
                wifi_params.MCS(8).Repetition = 1;
                wifi_params.MCS(8).M = 4; % pi/2-QPSK
                wifi_params.MCS(8).coding_type = 'LDPC';
                wifi_params.MCS(8).N_cbps = 2; % Number of coded bits per symbol
                wifi_params.MCS(8).N_cbpb = 2*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(8).DataRate = [1650.0, 1540.0, 1320.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 8 -------------------------------------------------------------------
                wifi_params.MCS(9).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(9).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(9).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(9).field_base = 8; % Base MCS field
                wifi_params.MCS(9).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(9).MCS = 'MCS8';
                wifi_params.MCS(9).CR = 5/8;
                wifi_params.MCS(9).Repetition = 1;
                wifi_params.MCS(9).M = 4; % pi/2-QPSK
                wifi_params.MCS(9).coding_type = 'LDPC';
                wifi_params.MCS(9).N_cbps = 2; % Number of coded bits per symbol
                wifi_params.MCS(9).N_cbpb = 2*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(9).DataRate = [2062.5, 1925.0, 1650.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 9 -------------------------------------------------------------------
                wifi_params.MCS(10).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(10).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(10).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(10).field_base = 9; % Base MCS field
                wifi_params.MCS(10).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(10).MCS = 'MCS9';
                wifi_params.MCS(10).CR = 3/4;
                wifi_params.MCS(10).Repetition = 1;
                wifi_params.MCS(10).M = 4; % pi/2-QPSK
                wifi_params.MCS(10).coding_type = 'LDPC';
                wifi_params.MCS(10).N_cbps = 2; % Number of coded bits per symbol
                wifi_params.MCS(10).N_cbpb = 2*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(10).DataRate = [2475.0, 2310.0, 1980.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 10 -------------------------------------------------------------------
                wifi_params.MCS(11).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(11).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(11).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(11).field_base = 10; % Base MCS field
                wifi_params.MCS(11).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(11).MCS = 'MCS10';
                wifi_params.MCS(11).CR = 13/16;
                wifi_params.MCS(11).Repetition = 1;
                wifi_params.MCS(11).M = 4; % pi/2-QPSK
                wifi_params.MCS(11).coding_type = 'LDPC';
                wifi_params.MCS(11).N_cbps = 2; % Number of coded bits per symbol
                wifi_params.MCS(11).N_cbpb = 2*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(11).DataRate = [2681.25, 2502.5, 2145.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 11 -------------------------------------------------------------------
                wifi_params.MCS(12).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(12).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(12).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(12).field_base = 11; % Base MCS field
                wifi_params.MCS(12).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(12).MCS = 'MCS11';
                wifi_params.MCS(12).CR = 7/8;
                wifi_params.MCS(12).Repetition = 1;
                wifi_params.MCS(12).M = 4; % pi/2-QPSK
                wifi_params.MCS(12).coding_type = 'LDPC';
                wifi_params.MCS(12).N_cbps = 2; % Number of coded bits per symbol
                wifi_params.MCS(12).N_cbpb = 2*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(12).DataRate = [2887.5, 2695.0, 2310.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                switch use8PSK
                    case false % pi/2-16-QAM used
                        % MCS 12 -------------------------------------------------------------------
                        wifi_params.MCS(13).phy_type = 'SC'; % stands for Single carrier
                        wifi_params.MCS(13).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                        wifi_params.MCS(13).N_CB = 1; % Possible number of bonded channels/defined in batch file
                        wifi_params.MCS(13).field_base = 12; % Base MCS field
                        wifi_params.MCS(13).field_extended = []; % Extended SC MCS indication field
                        wifi_params.MCS(13).MCS = 'MCS12';
                        wifi_params.MCS(13).CR = 1/2;
                        wifi_params.MCS(13).Repetition = 1;
                        wifi_params.MCS(13).M = 16; % pi/2-16QAM
                        wifi_params.MCS(13).coding_type = 'LDPC';
                        wifi_params.MCS(13).N_cbps = 4; % Number of coded bits per symbol
                        wifi_params.MCS(13).N_cbpb = 4*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                        wifi_params.MCS(13).DataRate = [3300.0, 3080.0, 2640.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                        % MCS 13 -------------------------------------------------------------------
                        wifi_params.MCS(14).phy_type = 'SC'; % stands for Single carrier
                        wifi_params.MCS(14).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                        wifi_params.MCS(14).N_CB = 1; % Possible number of bonded channels/defined in batch file
                        wifi_params.MCS(14).field_base = 13; % Base MCS field
                        wifi_params.MCS(14).field_extended = []; % Extended SC MCS indication field
                        wifi_params.MCS(14).MCS = 'MCS13';
                        wifi_params.MCS(14).CR = 5/8;
                        wifi_params.MCS(14).Repetition = 1;
                        wifi_params.MCS(14).M = 16; % pi/2-16QAM
                        wifi_params.MCS(14).coding_type = 'LDPC';
                        wifi_params.MCS(14).N_cbps = 4; % Number of coded bits per symbol
                        wifi_params.MCS(14).N_cbpb = 4*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                        wifi_params.MCS(14).DataRate = [4125.0, 3850.0, 3300.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                    case true  % pi/2-8-PSK used
                        % MCS 12 -------------------------------------------------------------------
                        wifi_params.MCS(13).phy_type = 'SC'; % stands for Single carrier
                        wifi_params.MCS(13).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                        wifi_params.MCS(13).N_CB = 1; % Possible number of bonded channels/defined in batch file
                        wifi_params.MCS(13).field_base = 12; % Base MCS field
                        wifi_params.MCS(13).field_extended = []; % Extended SC MCS indication field
                        wifi_params.MCS(13).MCS = 'MCS12';
                        wifi_params.MCS(13).CR = 2/3;
                        wifi_params.MCS(13).Repetition = 1;
                        wifi_params.MCS(13).M = 8; % pi/2-8PSK
                        wifi_params.MCS(13).coding_type = 'LDPC';
                        wifi_params.MCS(13).N_cbps = 3; % Number of coded bits per symbol
                        wifi_params.MCS(13).N_cbpb = 3*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                        wifi_params.MCS(13).DataRate = [3300.0, 3080.0, 2640.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                        % MCS 13 -------------------------------------------------------------------
                        wifi_params.MCS(14).phy_type = 'SC'; % stands for Single carrier
                        wifi_params.MCS(14).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                        wifi_params.MCS(14).N_CB = 1; % Possible number of bonded channels/defined in batch file
                        wifi_params.MCS(14).field_base = 13; % Base MCS field
                        wifi_params.MCS(14).field_extended = []; % Extended SC MCS indication field
                        wifi_params.MCS(14).MCS = 'MCS13';
                        wifi_params.MCS(14).CR = 5/6;
                        wifi_params.MCS(14).Repetition = 1;
                        wifi_params.MCS(14).M = 8; % pi/2-8PSK
                        wifi_params.MCS(14).coding_type = 'LDPC';
                        wifi_params.MCS(14).N_cbps = 3; % Number of coded bits per symbol
                        wifi_params.MCS(14).N_cbpb = 3*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                        wifi_params.MCS(14).DataRate = [4125.0, 3850.0, 3300.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                end
                % MCS 14 -------------------------------------------------------------------
                wifi_params.MCS(15).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(15).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(15).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(15).field_base = 14; % Base MCS field
                wifi_params.MCS(15).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(15).MCS = 'MCS14';
                wifi_params.MCS(15).CR = 3/4;
                wifi_params.MCS(15).Repetition = 1;
                wifi_params.MCS(15).M = 16; % pi/2-16QAM
                wifi_params.MCS(15).coding_type = 'LDPC';
                wifi_params.MCS(15).N_cbps = 4; % Number of coded bits per symbol
                wifi_params.MCS(15).N_cbpb = 4*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(15).DataRate = [4950.0, 4620.0, 3960.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 15 -------------------------------------------------------------------
                wifi_params.MCS(16).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(16).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(16).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(16).field_base = 15; % Base MCS field
                wifi_params.MCS(16).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(16).MCS = 'MCS15';
                wifi_params.MCS(16).CR = 13/16;
                wifi_params.MCS(16).Repetition = 1;
                wifi_params.MCS(16).M = 16; % pi/2-16QAM
                wifi_params.MCS(16).coding_type = 'LDPC';
                wifi_params.MCS(16).N_cbps = 4; % Number of coded bits per symbol
                wifi_params.MCS(16).N_cbpb = 4*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(16).DataRate = [5362.5, 5005.0, 4290.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 16 -------------------------------------------------------------------
                wifi_params.MCS(17).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(17).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(17).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(17).field_base = 16; % Base MCS field
                wifi_params.MCS(17).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(17).MCS = 'MCS16';
                wifi_params.MCS(17).CR = 7/8;
                wifi_params.MCS(17).Repetition = 1;
                wifi_params.MCS(17).M = 16; % pi/2-16QAM
                wifi_params.MCS(17).coding_type = 'LDPC';
                wifi_params.MCS(17).N_cbps = 4; % Number of coded bits per symbol
                wifi_params.MCS(17).N_cbpb = 4*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(17).DataRate = [5775.0, 5390.0, 4620.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 17 -------------------------------------------------------------------
                wifi_params.MCS(18).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(18).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(18).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(18).field_base = 17; % Base MCS field
                wifi_params.MCS(18).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(18).MCS = 'MCS17';
                wifi_params.MCS(18).CR = 1/2;
                wifi_params.MCS(18).Repetition = 1;
                wifi_params.MCS(18).M = 64; % pi/2-64QAM
                wifi_params.MCS(18).coding_type = 'LDPC';
                wifi_params.MCS(18).N_cbps = 6; % Number of coded bits per symbol
                wifi_params.MCS(18).N_cbpb = 6*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(18).DataRate = [4950.0, 4620.0, 3960.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 18 -------------------------------------------------------------------
                wifi_params.MCS(19).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(19).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(19).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(19).field_base = 18; % Base MCS field
                wifi_params.MCS(19).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(19).MCS = 'MCS18';
                wifi_params.MCS(19).CR = 5/8;
                wifi_params.MCS(19).Repetition = 1;
                wifi_params.MCS(19).M = 64; % pi/2-64QAM
                wifi_params.MCS(19).coding_type = 'LDPC';
                wifi_params.MCS(19).N_cbps = 6; % Number of coded bits per symbol
                wifi_params.MCS(19).N_cbpb = 6*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(19).DataRate = [6187.5, 5775.0, 4950.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 19 -------------------------------------------------------------------
                wifi_params.MCS(20).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(20).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(20).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(20).field_base = 19; % Base MCS field
                wifi_params.MCS(20).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(20).MCS = 'MCS19';
                wifi_params.MCS(20).CR = 3/4;
                wifi_params.MCS(20).Repetition = 1;
                wifi_params.MCS(20).M = 64; % pi/2-64QAM
                wifi_params.MCS(20).coding_type = 'LDPC';
                wifi_params.MCS(20).N_cbps = 6; % Number of coded bits per symbol
                wifi_params.MCS(20).N_cbpb = 6*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(20).DataRate = [7425.0, 6930.0, 5940.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 20 -------------------------------------------------------------------
                wifi_params.MCS(21).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(21).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(21).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(21).field_base = 20; % Base MCS field
                wifi_params.MCS(21).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(21).MCS = 'MCS20';
                wifi_params.MCS(21).CR = 13/16;
                wifi_params.MCS(21).Repetition = 1;
                wifi_params.MCS(21).M = 64; % pi/2-64QAM
                wifi_params.MCS(21).coding_type = 'LDPC';
                wifi_params.MCS(21).N_cbps = 6; % Number of coded bits per symbol
                wifi_params.MCS(21).N_cbpb = 6*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(21).DataRate = [8043.75, 7507.5, 6435.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                % MCS 21 -------------------------------------------------------------------
                wifi_params.MCS(22).phy_type = 'SC'; % stands for Single carrier 
                wifi_params.MCS(22).N_ss = 1; % Possible numbers of spatial stream/defined in batch file
                wifi_params.MCS(22).N_CB = 1; % Possible number of bonded channels/defined in batch file
                wifi_params.MCS(22).field_base = 21; % Base MCS field
                wifi_params.MCS(22).field_extended = []; % Extended SC MCS indication field
                wifi_params.MCS(22).MCS = 'MCS21';
                wifi_params.MCS(22).CR = 7/8;
                wifi_params.MCS(22).Repetition = 1;
                wifi_params.MCS(22).M = 64; % pi/2-64QAM
                wifi_params.MCS(22).coding_type = 'LDPC';
                wifi_params.MCS(22).N_cbps = 6; % Number of coded bits per symbol
                wifi_params.MCS(22).N_cbpb = 6*[480 448 384]; % Number of coded bits per symbol block - to be select according to GI length
                wifi_params.MCS(22).DataRate = [8662.5, 8085.0, 6930.0]*1e6; % Max. theoretical throughput, or data rate (for short, normal or long GI length)
                
         
end
%% Select true MCS indices
mcs_tol = 0.001;
ind_mcs_offset = 0;
switch wifi_standard
    case '802dot11ad'
        ad_MCS_vec_all = [0, 1:9, 9.1, 10:12, 12.1:0.1:12.6, 13:31].';
        if ischar(MCSvec)
            switch MCSvec
                case 'all'
                    MCSvec = ad_MCS_vec_all;
                case 'Ctrl'
                    MCSvec = 0;
                case 'all-SC'
                    MCSvec = [1:9, 9.1, 10:12, 12.1:12.6];
                case 'all-OFDM'
                    MCSvec = 13:24;
                case 'all-LPSC'
                    MCSvec = 25:31;
                otherwise
                    error('Unsupported option');
            end
        end
        
        ind_mcs_tmp1 = find(abs(ad_MCS_vec_all-MCSvec) < mcs_tol);
        ind_mcs_mod1 = mod(ind_mcs_tmp1, length(ad_MCS_vec_all));
        ind_mcs_mod1(ind_mcs_mod1 == 0) = length(ad_MCS_vec_all);
        i_mcs_vec_aux = (sort(ind_mcs_mod1 + ind_mcs_offset)).';
        
        if isempty(i_mcs_vec_aux)
            error('Chosen MCS undefined for IEEE 802.11ad SC mode (see load_wifi_params.m and mcs_definition.m)');
        end
    case '802dot11ay'
        ay_MCS_vec_all = [0,... % for Control
            1:21].'; % for SC
        if ischar(MCSvec)
            switch MCSvec
                case 'all'
                    MCSvec = ay_MCS_vec_all;
                case 'Ctrl'
                    MCSvec = 0;
                case 'all-SC'
                    MCSvec = 1:20;
%                 case 'all-OFDM'
%                     MCSvec = 13:24;
                otherwise
                    error('Unsupported option');
            end
        end
        
        ind_mcs_tmp1 = find(abs(ay_MCS_vec_all-MCSvec) < mcs_tol);
        ind_mcs_mod1 = mod(ind_mcs_tmp1, length(ay_MCS_vec_all));
        ind_mcs_mod1(ind_mcs_mod1 == 0) = length(ay_MCS_vec_all);
        i_mcs_vec_aux = (sort(ind_mcs_mod1 + ind_mcs_offset)).';
        
    case {'802dot11g', '802dot11ac', '802dot11af', '802dot11ah', '802dot11ax'}
        warning('Check in future (see load_wifi_params.m)');
    otherwise
        error('Undefined IEEE 802.11 standard (see load_wifi_params.m)');
end

