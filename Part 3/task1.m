% Define constants
L = 199;
l = 0:1:L;
K = 2048;
k_range = 0:1:K-1;

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
zp = bb_rece_data_172648_1474(pilot_index+1);

% Define the V matrix
V = exp(-1i*2*pi.*((pilot_index.*l)/K));

% Need to find the D matrix using the OFDM_PILOT
D = diag(OFDM_PILOT(pilot_index+1));

% Now calculate the time-domain channel
h_ls = (1/K_p)*V'*D'*zp;

% Now calculate the frequeny-domain channel estimate
freq_v = exp(-1i*2*pi.*(((k_range.').*l)/K));

H = freq_v*h_ls;