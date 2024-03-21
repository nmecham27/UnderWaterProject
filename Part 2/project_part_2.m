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

%pb_cor_length = 1:length(pb_cor);
figure(2);
hold on
title("Correlated passband data");
plot(lag, pb_cor); %292723
hold off
