clear variables;
clf;
K = 2048; %The total number of subcarriers
L = 200; %The number of zero-padded symbols
Fc = 24*10^3; %Carrier Frequency
B = 8 * 10^3; % Bandwidth
sampling_rate_256 = 256*10^3;
sampling_rate_192 = 192*10^3;
W = 24; %The number of ZP-OFDM symbols in one packet
lambda = 24; %The oversampling factor
v = 1.03; %Underwater velocity in m/s
c = 1500; %Speed of sound in m/s
T_tx = 8.2695; %Transmitted signal duration in seconds
figure_number = 1; %Variable to track figure number so we don't have to update everywhere

%Parameters used for the resampling of step 4
Lp = 24;
Ms = 256;
Ls = 192;

N = Lp*Ls-1;
h = Ls*fir1(N,1/Ms,kaiser(N+1,7.8562));
%*****

%Load the benchmark data
load("benchmark_rece_data_174623_1472.mat");
load("pilot_signal_for_synchronization.mat"); %Loads the data into OFDM_data_pre_old
load("itc_1007_compfilter.mat"); % stores vector in h_comp variable
load('ofdm_map.mat');
load("benchmark_parameter_174623_1472.mat");
load('test_rece_data_173048_1472.mat');

%Reassign the values loaded in so that we can more easily change to other
%sets of data

%Benchmark data.
%y_rec_data = rece_data_ofdm_bench;

%Real data
y_rec_data = rece_data_ofdm_test;

%Pilot signal
y_pilot_signal = OFDM_data_pre_old;

%Plot the original y data
rec_data_length = 1:length(y_rec_data);

figure(figure_number);
figure_number = figure_number + 1;
hold on
title("y original data");
plot(rec_data_length, y_rec_data);
hold off

%Filter out the noise from the passband signal
y_pb = bandpass(y_rec_data,[-4000+Fc,4000+Fc], sampling_rate_256);

%The mach number
a = v/c;

%Resample with the estimated mach number, a_hat.
%Plotting the passband data to try and estimate the
%T_rx.
x = 1:length(y_pb);
x = x.';
figure(figure_number);
figure_number = figure_number + 1;
hold on
title("Passband filtered data");
plot(x,y_pb);
hold off

%From the plot I see something like the following for the T_rx
%Hopefully this is only only place now that we have to do manual
%data input
%Uncomment for real data
sample_diff=2127860-12243;
T_rx = sample_diff/sampling_rate_256;
a_hat = (T_tx/T_rx)-1;

% Uncomment for benchmark data
%a_hat = a_est_174623;

%Resample the data with a_hat
y_pb_re = resample(y_pb,round(((1+a_hat)*10^5)),(10^5));

%Resample the data to match our transmitter sampling rate of 192 KHz
y_pb_re_192 = upfirdn(y_pb_re, h, Ls, Ms);

y_pb_re_192_length = 1:length(y_pb_re_192);

%Plot the resampled data
figure(figure_number);
figure_number = figure_number + 1;
hold on
title("Resampled passband data (192Hz)");
plot(y_pb_re_192_length, y_pb_re_192);
hold off

%Plot the pilot signal
y_pilot_signal_length = 1:length(y_pilot_signal);
figure(figure_number);
figure_number = figure_number + 1;
hold on
title("Pilot signal");
plot(y_pilot_signal_length, y_pilot_signal);
hold off

%Perform the correlation between the filtered passband data and the
%pilot signal
[pb_cor,lag] = xcorr(y_pb_re_192, y_pilot_signal);
pb_cor_length = 1:length(pb_cor);
%Plot the correlation results
figure(figure_number);
figure_number = figure_number + 1;
hold on
title("Correlated passband data");
% No matter how you plot, peak is 1418947 samples before end of data
plot(lag, pb_cor); %292723
%plot(pb_cor_length, pb_cor); % Peak is at 2004390
hold off

%COULD POSSIBLY REMOVE THIS
% Want last 1418947 samples
% y_pb_re_192 is 1711669 samples long
%n0 = length(y_pb_re_192)-1418947;
%y_pb_re_192 = y_pb_re_192(n0:length(y_pb_re_192));
%**************************

%We can use the max value & index from the correlation results to find the
%index
[max_correlation, max_correlation_index] = max(pb_cor, [], "all", "linear");
n0 = lag(max_correlation_index);
y_pb_re_192 = y_pb_re_192(n0:length(y_pb_re_192));

%Plot the newly sized passband data
y_pb_re_192_length = 1:length(y_pb_re_192);
figure(figure_number);
figure_number = figure_number + 1;
hold on;
title("Passband data (192Hz) after resizing due to correlation results");
plot(y_pb_re_192_length, y_pb_re_192);
hold off;

% Step 6 (is this step necessary)
y_pb_re_192_conv = conv(y_pb_re_192, h_comp); % adds 50 sample delay
figure(figure_number);
figure_number = figure_number + 1;
hold on
title("y_pb_re_192_conv");
plot(1:length(y_pb_re_192_conv), y_pb_re_192_conv);
hold off
% Need to remove 50 sample delay
y_pb_re_192_conv = y_pb_re_192_conv(51:length(y_pb_re_192_conv)); % Not sure if this is correct way to remove delay

% Step 7
% Convert from passband to baseband
yBB_I = 1:length(y_pb_re_192_conv);
yBB_Q = 1:length(y_pb_re_192_conv);
yBB_I = yBB_I.';
yBB_Q = yBB_Q.';
yBB = 1:length(y_pb_re_192_conv);
yBB = yBB.';
ts_192 = 1/sampling_rate_192;
for n = 1:length(yBB_I)
    yBB_I(n) = y_pb_re_192_conv(n)*2*cos(2*pi*Fc*n*ts_192);
    yBB_Q(n) = -1*y_pb_re_192_conv(n)*2*sin(2*pi*Fc*n*ts_192);
    yBB(n) = yBB_I(n) + j*yBB_Q(n);
end

% Step 8
% Raised cosine filter (arguments come from Part 1 of project)
beta = .125;
delay = 100;
filter = rcosine(1, lambda, 'sqrt', beta, delay);
filter = filter.';
yBB_filtered = conv(yBB, filter);
% Remove delay from filtered data
yBB_filtered = yBB_filtered((lambda*delay)+1:length(yBB_filtered));
yBB_filtered = yBB_filtered.';

% Step 9
% Grid search to find the starting index of the 1st OFDM block and the
% carrier frequency offset of the 1st OFDM block

n_index = 0:1:(K+L)*lambda-1;
down_sample_index = 0:1:K+L-1;

%***************************
%FFT matrix
m = [1:K].';
j = [1:K+L];
fft_matrix = exp(-1i*((2*pi*(m-1)*(j-1))/(K)));
%***************************


p_null = zeros(length(2200:1:2400), length(-2:.1:2));
%Loop through the sample index
%Note: This loop takes a LONG time. Like 10 minutes to run
p_null_n_index = 1;
p_null_eps_index = 1;
for n = 2200:1:2400
    % Loop through the frequency offsets
    for eps = -2:.1:2
        yBB_cfo_comp = zeros(1,(K+L)*lambda);

        yBB_cfo_comp(n_index+1) = yBB_filtered(n_index+n+1).*exp(-1i*2*pi*eps*(n_index+n+1)*ts_192);
        
        %Down sampling
        yBB_down_sampled = zeros(1, K+L);
        yBB_down_sampled(down_sample_index+1) = yBB_cfo_comp((down_sample_index*lambda)+1);
        yBB_down_sampled = yBB_down_sampled.';

        %Obtaining frequency data
        %This is the biggest problem loop for why everything takes so long
        z_freq_data = zeros(K,1);
        z_freq_data = fft_matrix*yBB_down_sampled;


        %Calculate the power over null subcarriers
        %Use the ofdm_map data provided to find the index of the null sub
        %carriers.
        p_null(p_null_n_index, p_null_eps_index) = sum(abs(z_freq_data(ofdm_map==0)).^2);
        p_null_eps_index = p_null_eps_index + 1;
    end
    p_null_eps_index = 1;
    p_null_n_index = p_null_n_index + 1;
end

%Benchmark save
%save("EndOfPart9.mat", "p_null", "z_freq_data", "yBB_filtered", "K", "lambda", "L", "a_hat");
save("EndOfPart9Real.mat", "p_null", "z_freq_data", "yBB_filtered", "K", "lambda", "L", "a_hat");