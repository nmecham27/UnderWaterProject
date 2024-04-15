% Define constants
L = 199;
l = 0:1:L;
K = 2048;
k_range = 0:1:K-1;
W = 21;

%% Part a
% Load the test data
load("OFDM_PILOT.mat");
load("ofdm_map.mat");
load("benchmark_parameter_172648_1.mat");
load("benchmark_Zw_172648_1.mat");

% Obtain the indices of the pilot symbols
% NOTE: The indices found in the vector are actual_index-1
pilot_index = find(ofdm_map == 1)-1;
K_p = length(pilot_index);

% The Zw results but only the data corresponding to
% pilot subcarriers
zp1 = bb_rece_data_172648_1474(pilot_index+1,:);

% Define the V matrix
V = exp(-1i*2*pi.*((pilot_index.*l)/K));

% Need to find the D matrix using the OFDM_PILOT
D = diag(OFDM_PILOT(pilot_index+1));

% Now calculate the time-domain channel
h_ls_one = (1/K_p)*V'*D'*zp1;

% Now calculate the frequeny-domain channel estimate
freq_v = exp(-1i*2*pi.*(((k_range.').*l)/K));

H1 = freq_v*h_ls_one;

%% Part b
% Find the indices of the null subcarriers
null_index = find(ofdm_map == 0)-1;
K_n = length(null_index);

% Using the null subcarriers calculate the noise variance
noise_variance = (1/K_n)*sum(abs(bb_rece_data_172648_1474(null_index+1,:)).^2);

% Find the indices of the data subcarriers
data_index = find(ofdm_map == 2)-1;
K_d = length(data_index);

% The Zw results but only the data corresponding to
% data subcarriers
zd = bb_rece_data_172648_1474(data_index+1,:);

% Get the channel estimate corresponding to the data
% subcarriers
H_d = H1(data_index+1,:);

% Calculate the LR
% The format for the LR matrix is each row corresponds to a
% differnt OFDM symbol and the colums are L_b1 and L_b1 for
% each subcarrier
LR = zeros(2*K_d, W);
for ofdm_symbol = 1:W
    for subcarrier = 1:K_d
%         x_1 = -1*(norm(zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol)*(1/sqrt(2)+1i*1/sqrt(2)))^2)/noise_variance(ofdm_symbol);
%         x_2 = -1*(norm(zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol)*(-1/sqrt(2)+1i*1/sqrt(2)))^2)/noise_variance(ofdm_symbol);
%         x_3 = -1*(norm(zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol)*(1/sqrt(2)-1i*1/sqrt(2)))^2)/noise_variance(ofdm_symbol);
%         x_4 = -1*(norm(zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol)*(-1/sqrt(2)-1i*1/sqrt(2)))^2)/noise_variance(ofdm_symbol);
% 
%         L_b1 = (max(x_1,x_3)*log(1+exp(-1*abs(x_3-x_1)))) - (max(x_2,x_4)*log(1+exp(-1*abs(x_4-x_2))));
%         L_b2 = (max(x_1,x_2)*log(1+exp(-1*abs(x_2-x_1)))) - (max(x_3,x_4)*log(1+exp(-1*abs(x_4-x_3))));
% 
%         LR(ofdm_symbol, (subcarrier*2)-1) = L_b1;
%         LR(ofdm_symbol, subcarrier*2) = L_b2;

        x1 = 1/sqrt(2) + i*1/sqrt(2);
        x2 = -1/sqrt(2) + i*1/sqrt(2);
        x3 = 1/sqrt(2) - i*1/sqrt(2);
        x4 = -1/sqrt(2) - i*1/sqrt(2);
        
        A = norm(-1*zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol).*x1)^2/noise_variance(ofdm_symbol);
        B = norm(-1*zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol).*x3)^2/noise_variance(ofdm_symbol);
        C = norm(-1*zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol).*x2)^2/noise_variance(ofdm_symbol);
        D = norm(-1*zd(subcarrier,ofdm_symbol)-H_d(subcarrier,ofdm_symbol).*x4)^2/noise_variance(ofdm_symbol);
        
        L_b1 = max(A, B)+log(1+exp(-1*abs(B-A)))-(max(C,D)+log(1+exp(-1*abs(D-C))));
        L_b2 = max(A, C)+log(1+exp(-1*abs(C-A)))-(max(B,D)+log(1+exp(-1*abs(D-B))));
        
        LR((subcarrier*2)-1, ofdm_symbol) = L_b1;
        LR(subcarrier*2, ofdm_symbol) = L_b2;
    end
end

%% Part c
% Load interleaver vector and codewords
load("INTRLVR.mat");
load("CODE.mat");

NAME = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt';
[address, LDPC_INFOLEN] = ldpc_mex_initial_CAPI([1420, 2840, 2], NAME);
est_code = zeros(2*K_d, W);
bec = zeros(1, W);
LR_in_de = zeros(length(LR), 1);
APP_code = zeros(2*K_d, W);
wec = 0;
for ofdm_sym = 1:W
    LR_in_de(INTRLVR) = LR(:,ofdm_sym);
    APP_code(:,ofdm_sym) = ldpcDecoder_CAPI(address, LR_in_de);
    est_code(:,ofdm_sym) = (APP_code(:,ofdm_sym) < 0);
    count = sum(abs(est_code(:,ofdm_sym)-CODE(:,ofdm_sym)));
    bec(ofdm_sym) = count;
    if count ~= 0
        wec = wec + 1;
    end
end

ber = bec/W/length(APP_code);
bler = wec/W;


%% Task 2
%% Part a
% Hydrophone 2
load("benchmark_Zw_172648_2.mat");

zp2 = bb_rece_data_172648_1475(pilot_index+1,:);

% Need to redefine D matrix, D variable was redefined when calculating LLR
D = diag(OFDM_PILOT(pilot_index+1));

% Now calculate the time-domain channel
h_ls_two = (1/K_p)*V'*D'*zp2;
% Freq domain channel
H2 = freq_v*h_ls_two;

% Hydrophone 3
load("benchmark_Zw_172648_3.mat");

zp3 = bb_rece_data_172648_1476(pilot_index+1,:);
% Now calculate the time-domain channel
h_ls_three = (1/K_p)*V'*D'*zp3;
% Freq domain channel
H3 = freq_v*h_ls_three;

%% Part b

