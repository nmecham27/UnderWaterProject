clear all;
clf;
K = 2048; %The total number of subcarriers
L = 200; %The number of zero-padded symbols
Fc = 24*10^3; %Carrier Frequency
B = 8 * 10^3; % Bandwidth
sampling_rate_high = 256*10^3;
sampling_rate_low = 192*10^3;
W = 24; %The number of ZP-OFDM symbols in one packet
lambda = 24; %The oversampling factor
v = 1.03; %Underwater velocity in m/s
c = 1500; %Speed of sound in m/s
T_tx = 8.2695; %Transmitted signal duration in seconds

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

%Filter out the noise from the passband signal
y_pb = bandpass(rece_data_ofdm_bench,[-1000+Fc,8000+Fc], sampling_rate_high);

%The mach number
a = v/c;

%Resample with the estimated mach number, a_hat.
%Plotting the passband data to try and estimate the
%T_rx.
x = 1:length(y_pb);
x = x.';
figure(1);
hold on
title("Passband filtered data");
plot(x,y_pb);
hold off

%From the plot I see something like the following for the T_rx
sample_diff=2121170-4269;
T_rx = sample_diff/sampling_rate_high;

a_hat = (T_tx/T_rx)-1;

%Resample the data with a_hat
y_pb_re = resample(y_pb,round(((1+a_hat)*10^5)),(10^5));

%Resample the data to match our transmitter sampling rate of 192 KHz
y_pb_re_192 = upfirdn(y_pb_re, h, Ls, Ms);

[pb_cor,lag] = xcorr(y_pb_re_192, OFDM_data_pre_old);

pb_cor_length = 1:length(pb_cor);
figure(2);
hold on
title("Correlated passband data");
% No matter how you plot, peak is 1418947 samples before end of data
plot(lag, pb_cor); %292723
%plot(pb_cor_length, pb_cor); % Peak is at 2004390
hold off

% Want last 1418947 samples
% y_pb_re_192 is 1711669 samples long
n0 = length(y_pb_re_192)-1418947;
y_pb_re_192 = y_pb_re_192(n0:length(y_pb_re_192));

% Step 6 (is this step necessary)
load("itc_1007_compfilter.mat"); % stores vector in h_comp variable
y_pb_re_192_conv = conv(y_pb_re_192, h_comp); % adds 50 sample delay
figure(3);
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
ts_low = 1/sampling_rate_low;
for n = 1:length(yBB_I)
    yBB_I(n) = y_pb_re_192_conv(n)*2*cos(2*pi*Fc*n*ts_low);
    yBB_Q(n) = -1*y_pb_re_192_conv(n)*2*sin(2*pi*Fc*n*ts_low);
    yBB(n) = yBB_I(n) + j*yBB_Q(n);
end

% Step 8
% Raised cosine filter (arguments come from Part 1 of project)
lambda = 24;
beta = .125;
delay = 100;
filter = rcosine(1, lambda, 'sqrt', beta, delay);
filter = filter.';
yBB_filtered = conv(yBB, filter);
% Remove delay from filtered data
yBB_filtered = yBB_filtered((lambda*delay*2)+1:length(yBB_filtered));


