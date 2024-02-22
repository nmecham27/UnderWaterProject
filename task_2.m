clear all;
clf;
M = 2048; % # of subcarriers (2048)
L = 201; % # of channel taps
N = 1; % # of OFDM symbols (16)
h = [0.227 0.46 0.688 0.46 0.227];
lambda = 24;
delay = 100;
beta = .125;
fs = 192e3;
B = 8000; % bandwidth is 8000 Hz
Ts = 1/B;
delta_f = B/M;
T = 1/delta_f;
Tg = .025; 
fc = 24e3;
ts = Ts/lambda;


original_symbols = [];
x = [];
for i = 1:N
    % Generate M QPSK symbols (Step 1)
    data = randi([0 1], 2*M, 1);
    data_freq = nrSymbolModulate(data, 'QPSK');
    
    original_symbols = [original_symbols; data_freq];

    % Convert to time domain (Step 2)
    data_time = ifft(data_freq);
    
    % Zero padding (Step 3)
    x = [x; zeros(L-1, 1); data_time];
end
% Ending 0's (Step 3)
x = [x; zeros(L-1, 1)];

% Upsample (step 4)
x = upsample(x, lambda);

% Pre-filter Plots
f = (fs*(-length(x)/2:length(x)/2 - 1)/length(x));
figure(1)
plot(f, abs(x));
hold on
title("Pre filtered samples");
hold off
figure(2)
spectrogram(abs(x),[],[],[],fs,'yaxis');
%spectrogram(x,[],[],[],fs,'yaxis');
hold on
title("Pre filtered samples spectrogram");
hold off

% Raised cosine filter (Step 5)
filter = rcosine(1, lambda, 'sqrt', beta, delay);
filter = filter.';
x_filtered = conv(x, filter);

% Post-filter Plots
f_filtered = (fs*(-length(x_filtered)/2:length(x_filtered)/2 - 1)/length(x_filtered));
figure(3)
plot(f_filtered, abs(x_filtered));
hold on
title("Filtered samples");
hold off
figure(4)
spectrogram(abs(f_filtered),100,0,100,fs,'yaxis');
%spectrogram(f_filtered,100,0,100,fs,'yaxis');
hold on
title("Filtered samples spectrogram");
hold off

% Convert to pass band (Step 6)
x_pass_band = zeros(length(x_filtered), 1);
for i = 0:length(x_filtered)-1
    x_pass_band(i+1) = real(x_filtered(i+1)*exp(j*2*pi*fc*i*ts));
end

% LFM (Step 7)
f0 = fc - 4e3;
f1 = fc + 4e3;
t1 = .05;
b = (f1-f0)/t1;
t = linspace(0, t1, t1*fs);
f_t = f0+b.*t;
chirp = cos(2*pi.*f_t.*t).';
z = zeros(4*length(t), 1);
LFM = [];
for i = 1:4
   LFM = [LFM; chirp; z]; 
end
% LFM with OFDM symbols
lfm_symbols = [LFM; x_pass_band; LFM];

figure(5)
spectrogram(x_pass_band, 4096, (3/4)*4096, 4096, fs, 'yaxis');
hold on
title("OFDM passband signal")
hold off

figure(6)
spectrogram(LFM, 4096, (3/4)*4096, 4096, fs, 'yaxis');
hold on
title("Chirps")
hold off

figure(7)
spectrogram(lfm_symbols, 100, 80, 100, fs, 'yaxis');
hold on
title("OFDM symbols with chirp applied")
hold off

% Validation (using only passband data)
% Convert to baseband
lfm_symbols_base_band_I = zeros(length(x_pass_band), 1);
lfm_symbols_base_band_Q = zeros(length(x_pass_band), 1);
for n = 0: length(x_pass_band)-1
    lfm_symbols_base_band_I(n+1) = x_pass_band(n+1)*2*cos(2*pi*fc*n*ts); 
    lfm_symbols_base_band_Q(n+1) = -1*x_pass_band(n+1)*2*sin(2*pi*fc*n*ts); 
end
% Combine I and Q data for complex values
lfm_symbols_base_band = lfm_symbols_base_band_I + j*lfm_symbols_base_band_Q;

% Pass through filter
x_ovf = conv(lfm_symbols_base_band, filter);

% Remove delay due to filter
x_ovf = x_ovf((lambda*delay*2)+1:length(x_ovf));
% Downsample
x_ovf_down_sampled = downsample(x_ovf, lambda);
% Remove leading 0's due to zero padding
x_ovf_down_sampled = x_ovf_down_sampled(L:length(x_ovf_down_sampled));
% Remove trailing 0's due to zero padding
x_ovf_down_sampled = x_ovf_down_sampled(1:M);
% Go back to frequency domain
result = fft(x_ovf_down_sampled);
