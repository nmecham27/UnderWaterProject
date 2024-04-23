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

load("FinalData1.mat");
bb_rece_data_172648_1474 = z_results;

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
noise_variance_one = (1/K_n)*sum(abs(bb_rece_data_172648_1474(null_index+1,:)).^2);

% Find the indices of the data subcarriers
data_index = find(ofdm_map == 2)-1;
K_d = length(data_index);

% The Zw results but only the data corresponding to
% data subcarriers
zd_one = bb_rece_data_172648_1474(data_index+1,:);

% Get the channel estimate corresponding to the data
% subcarriers
H_d_one = H1(data_index+1,:);

% Calculate the LR
% The format for the LR matrix is each row corresponds to a
% differnt OFDM symbol and the colums are L_b1 and L_b1 for
% each subcarrier
LR = zeros(2*K_d, W);
x1 = 1/sqrt(2) + i*1/sqrt(2);
x2 = -1/sqrt(2) + i*1/sqrt(2);
x3 = 1/sqrt(2) - i*1/sqrt(2);
x4 = -1/sqrt(2) - i*1/sqrt(2);
for ofdm_symbol = 1:W
    for subcarrier = 1:K_d
        
        A = norm(-1*zd_one(subcarrier,ofdm_symbol)-H_d_one(subcarrier,ofdm_symbol).*x1)^2/noise_variance_one(ofdm_symbol);
        B = norm(-1*zd_one(subcarrier,ofdm_symbol)-H_d_one(subcarrier,ofdm_symbol).*x3)^2/noise_variance_one(ofdm_symbol);
        C = norm(-1*zd_one(subcarrier,ofdm_symbol)-H_d_one(subcarrier,ofdm_symbol).*x2)^2/noise_variance_one(ofdm_symbol);
        D = norm(-1*zd_one(subcarrier,ofdm_symbol)-H_d_one(subcarrier,ofdm_symbol).*x4)^2/noise_variance_one(ofdm_symbol);
        
        L_b1 = max(A, B)+log(1+exp(-1*abs(B-A)))-(max(C,D)+log(1+exp(-1*abs(D-C))));
        L_b2 = max(A, C)+log(1+exp(-1*abs(C-A)))-(max(B,D)+log(1+exp(-1*abs(D-B))));
        
        LR((subcarrier*2)-1, ofdm_symbol) = L_b1;
        LR(subcarrier*2, ofdm_symbol) = L_b2;
    end
end

% Scale down values for decoder
for ofdm_sym = 1:W
    if max(abs(LR(:,ofdm_sym)))>100
        max_val=max(LR(:,ofdm_sym));
        min_val=min(LR(:,ofdm_sym));
        if abs(max_val)>abs(min_val)
            LR(:,ofdm_sym)=LR(:,ofdm_sym)/abs(max_val)*100;
        else
            LR(:,ofdm_sym)=LR(:,ofdm_sym)/abs(min_val)*100;
        end
    end
end

%% Part c
% Load interleaver vector and codewords
load("INTRLVR.mat");
load("CODE.mat");

NAME = './5G_LDPC_M10_N20_Z142_Q2_nonVer.txt';
[address, LDPC_INFOLEN] = ldpc_mex_initial_CAPI([1420, 2840, 2], NAME);
est_code_one = zeros(2*K_d, W);
bec_one = zeros(1, W);
LR_in_de = zeros(length(LR), 1);
APP_code = zeros(2*K_d, W);
wec_one = 0;
for ofdm_sym = 1:W
    LR_in_de(INTRLVR) = LR(:,ofdm_sym);
    APP_code(:,ofdm_sym) = ldpcDecoder_CAPI(address, LR_in_de);
    est_code_one(:,ofdm_sym) = (APP_code(:,ofdm_sym) < 0);
    count = sum(abs(est_code_one(:,ofdm_sym)-CODE(:,ofdm_sym)));
    bec_one(ofdm_sym) = count;
    if count ~= 0 && ofdm_sym > 1
        wec_one = wec_one + 1;
    end
end

%ber_one = bec_one/W/length(APP_code); % Shouldnt divide by W here?
% Only calculate ber and bler for OFDM symbols 2-21
ber_one_total = sum(bec_one(:, 2:W))/(W-1)/length(APP_code);
bler_one = wec_one/(W-1);


%% Task 2
%% Part a
% Hydrophone 2
load("FinalData2.mat");
bb_rece_data_172648_1475 = z_results;

zp2 = bb_rece_data_172648_1475(pilot_index+1,:);

% Need to redefine D matrix, D variable was redefined when calculating LLR
D = diag(OFDM_PILOT(pilot_index+1));

% Now calculate the time-domain channel
h_ls_two = (1/K_p)*V'*D'*zp2;
% Freq domain channel
H2 = freq_v*h_ls_two;

% Hydrophone 3
load("FinalData3.mat");
bb_rece_data_172648_1476 = z_results;

zp3 = bb_rece_data_172648_1476(pilot_index+1,:);
% Now calculate the time-domain channel
h_ls_three = (1/K_p)*V'*D'*zp3;
% Freq domain channel
H3 = freq_v*h_ls_three;

%% Part b
% Get the channel estimate corresponding to the data
% subcarriers of second hydrophone
H_d_two = H2(data_index+1,:);

% Calculate norm(Hk) using only subcarriers data was sent
norm_Hk_2_hydro = zeros(K_d, W);
for ofdm_sym = 1: W
    for subc = 1: K_d
        norm_Hk_2_hydro(subc, ofdm_sym) =  sqrt((H_d_one(subc, ofdm_sym)*H_d_one(subc, ofdm_sym)')+(H_d_two(subc, ofdm_sym)*H_d_two(subc, ofdm_sym)')); 
    end
end

% Received data on hydrophone 2
zd_two = bb_rece_data_172648_1475(data_index+1,:);
% Variance on hydrophone 2
noise_variance_two = (1/K_n)*sum(abs(bb_rece_data_172648_1475(null_index+1,:)).^2);

% Calculate variance when combining hydrophones 1 and 2
noise_variance_2_hydro = zeros(K_d, W);
for ofdm_sym = 1:W
    for subc = 1: K_d
        noise_variance_2_hydro(subc, ofdm_sym) =  (1/(norm_Hk_2_hydro(subc, ofdm_sym)^2))*(noise_variance_one(ofdm_sym)*(H_d_one(subc, ofdm_sym)*H_d_one(subc, ofdm_sym)')+noise_variance_two(ofdm_sym)*(H_d_two(subc, ofdm_sym)*H_d_two(subc, ofdm_sym)'));
    end
end

% Calculate qk (scalar)
qk_2_hydro = zeros(K_d, W);
for ofdm_sym = 1: W
   for subc = 1:K_d
      qk_2_hydro(subc, ofdm_sym) = ([H_d_one(subc, ofdm_sym);H_d_two(subc, ofdm_sym)]'*[zd_one(subc, ofdm_sym);zd_two(subc, ofdm_sym)])/norm_Hk_2_hydro(subc, ofdm_sym); 
   end
end

% Calculate LLR
LR_2_hydro = zeros(2*K_d, W);
for ofdm_sym = 1:W
    for subc = 1:K_d
        A = (-1*abs(qk_2_hydro(subc, ofdm_sym)-norm_Hk_2_hydro(subc, ofdm_sym)*x1)^2)/noise_variance_2_hydro(subc, ofdm_sym);
        B = (-1*abs(qk_2_hydro(subc, ofdm_sym)-norm_Hk_2_hydro(subc, ofdm_sym)*x3)^2)/noise_variance_2_hydro(subc, ofdm_sym);
        C = (-1*abs(qk_2_hydro(subc, ofdm_sym)-norm_Hk_2_hydro(subc, ofdm_sym)*x2)^2)/noise_variance_2_hydro(subc, ofdm_sym);
        D = (-1*abs(qk_2_hydro(subc, ofdm_sym)-norm_Hk_2_hydro(subc, ofdm_sym)*x4)^2)/noise_variance_2_hydro(subc, ofdm_sym);
        
        L_b1 = max(A, B)+log(1+exp(-1*abs(B-A)))-(max(C,D)+log(1+exp(-1*abs(D-C))));
        L_b2 = max(A, C)+log(1+exp(-1*abs(C-A)))-(max(B,D)+log(1+exp(-1*abs(D-B))));
        
        LR_2_hydro((subc*2)-1, ofdm_sym) = L_b1;
        LR_2_hydro(subc*2, ofdm_sym) = L_b2;
    end
end

% Scale down values for decoder
for ofdm_sym = 1:W
    if max(abs(LR_2_hydro(:,ofdm_sym)))>100
        max_val=max(LR_2_hydro(:,ofdm_sym));
        min_val=min(LR_2_hydro(:,ofdm_sym));
        if abs(max_val)>abs(min_val)
            LR_2_hydro(:,ofdm_sym)=LR_2_hydro(:,ofdm_sym)/abs(max_val)*100;
        else
            LR_2_hydro(:,ofdm_sym)=LR_2_hydro(:,ofdm_sym)/abs(min_val)*100;
        end
    end
end

%% Part c
% Calculate ber and bler with 2 hydrophones
est_code_two = zeros(2*K_d, W);
bec_two = zeros(1, W);
LR_in_de = zeros(length(LR_2_hydro), 1);
APP_code = zeros(2*K_d, W);
wec_two = 0;
for ofdm_sym = 1:W
    LR_in_de(INTRLVR) = LR_2_hydro(:,ofdm_sym);
    APP_code(:,ofdm_sym) = ldpcDecoder_CAPI(address, LR_in_de);
    est_code_two(:,ofdm_sym) = (APP_code(:,ofdm_sym) < 0);
    count = sum(abs(est_code_two(:,ofdm_sym)-CODE(:,ofdm_sym)));
    bec_two(ofdm_sym) = count;
    if count ~= 0 && ofdm_sym > 1
        wec_two = wec_two + 1;
    end
end

%ber_two = bec_two/W/length(APP_code); % Shouldnt divide by W here?
ber_two_total = sum(bec_two(:, 2:W))/(W-1)/length(APP_code);
bler_two = wec_two/(W-1);

%% Part d
% Get the channel estimate corresponding to the data
% subcarriers of third hydrophone
H_d_three = H3(data_index+1,:);

% Calculate norm(Hk) using only subcarriers data was sent
norm_Hk_3_hydro = zeros(K_d, W);
for ofdm_sym = 1: W
    for subc = 1: K_d
        norm_Hk_3_hydro(subc, ofdm_sym) =  sqrt((H_d_one(subc, ofdm_sym)*H_d_one(subc, ofdm_sym)')+(H_d_two(subc, ofdm_sym)*H_d_two(subc, ofdm_sym)')+(H_d_three(subc, ofdm_sym)*H_d_three(subc, ofdm_sym)'));
    end
end


% Received data on hydrophone 3
zd_three = bb_rece_data_172648_1476(data_index+1,:);
% Variance on hydrophone 3
noise_variance_three = (1/K_n)*sum(abs(bb_rece_data_172648_1476(null_index+1,:)).^2);

% Calculate variance when combining hydrophones 1, 2, and 3
noise_variance_3_hydro = zeros(K_d, W);
for ofdm_sym = 1:W
    for subc = 1: K_d
        noise_variance_3_hydro(subc, ofdm_sym) =  (1/(norm_Hk_2_hydro(subc, ofdm_sym)^2))*(noise_variance_one(ofdm_sym)*(H_d_one(subc, ofdm_sym)*H_d_one(subc, ofdm_sym)')+noise_variance_two(ofdm_sym)*(H_d_two(subc, ofdm_sym)*H_d_two(subc, ofdm_sym)')+noise_variance_three(ofdm_sym)*(H_d_three(subc, ofdm_sym)*H_d_three(subc, ofdm_sym)'));
    end
end

% Calculate qk (scalar)
qk_3_hydro = zeros(K_d, W);
for ofdm_sym = 1: W
   for subc = 1:K_d
      qk_3_hydro(subc, ofdm_sym) = ([H_d_one(subc, ofdm_sym);H_d_two(subc, ofdm_sym);H_d_three(subc, ofdm_sym)]'*[zd_one(subc, ofdm_sym);zd_two(subc, ofdm_sym);zd_three(subc, ofdm_sym)])/norm_Hk_3_hydro(subc, ofdm_sym); 
   end
end

% Calculate LLR
LR_3_hydro = zeros(2*K_d, W);
for ofdm_sym = 1:W
    for subc = 1:K_d
        A = (-1*abs(qk_3_hydro(subc, ofdm_sym)-norm_Hk_3_hydro(subc, ofdm_sym)*x1)^2)/noise_variance_3_hydro(subc, ofdm_sym);
        B = (-1*abs(qk_3_hydro(subc, ofdm_sym)-norm_Hk_3_hydro(subc, ofdm_sym)*x3)^2)/noise_variance_3_hydro(subc, ofdm_sym);
        C = (-1*abs(qk_3_hydro(subc, ofdm_sym)-norm_Hk_3_hydro(subc, ofdm_sym)*x2)^2)/noise_variance_3_hydro(subc, ofdm_sym);
        D = (-1*abs(qk_3_hydro(subc, ofdm_sym)-norm_Hk_3_hydro(subc, ofdm_sym)*x4)^2)/noise_variance_3_hydro(subc, ofdm_sym);
        
        L_b1 = max(A, B)+log(1+exp(-1*abs(B-A)))-(max(C,D)+log(1+exp(-1*abs(D-C))));
        L_b2 = max(A, C)+log(1+exp(-1*abs(C-A)))-(max(B,D)+log(1+exp(-1*abs(D-B))));
        
        LR_3_hydro((subc*2)-1, ofdm_sym) = L_b1;
        LR_3_hydro(subc*2, ofdm_sym) = L_b2;
    end
end

% Scale down values for decoder
for ofdm_sym = 1:W
    if max(abs(LR_3_hydro(:,ofdm_sym)))>100
        max_val=max(LR_3_hydro(:,ofdm_sym));
        min_val=min(LR_3_hydro(:,ofdm_sym));
        if abs(max_val)>abs(min_val)
            LR_3_hydro(:,ofdm_sym)=LR_3_hydro(:,ofdm_sym)/abs(max_val)*100;
        else
            LR_3_hydro(:,ofdm_sym)=LR_3_hydro(:,ofdm_sym)/abs(min_val)*100;
        end
    end
end

% Calculate ber and bler with 3 hydrophones
est_code_three = zeros(2*K_d, W);
bec_three = zeros(1, W);
LR_in_de = zeros(length(LR_3_hydro), 1);
APP_code = zeros(2*K_d, W);
wec_three = 0;
for ofdm_sym = 1:W
    LR_in_de(INTRLVR) = LR_3_hydro(:,ofdm_sym);
    APP_code(:,ofdm_sym) = ldpcDecoder_CAPI(address, LR_in_de);
    est_code_three(:,ofdm_sym) = (APP_code(:,ofdm_sym) < 0);
    count = sum(abs(est_code_three(:,ofdm_sym)-CODE(:,ofdm_sym)));
    bec_three(ofdm_sym) = count;
    if count ~= 0 && ofdm_sym > 1
        wec_three = wec_three + 1;
    end
end

%ber_three = bec_three/W/length(APP_code); % Shouldnt divide by W here?
ber_three_total = sum(bec_three(:, 2:W))/(W-1)/length(APP_code);
bler_three = wec_three/(W-1);

save("part3_results.mat", "ber_one_total", "ber_two_total", "ber_three_total", "bler_one", "bler_two", "bler_three")