M = 8; % # of subcarriers (2048)
L = 5; % # of channel taps
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
    
    
    % Generate M QPSK symbols
    data = randi([0 1], 2*M, 1);
    data_freq = nrSymbolModulate(data, 'QPSK');
    
    original_symbols = [original_symbols; data_freq];

    % Convert to time domain
    data_time = ifft(data_freq);
    
    x = [x; zeros(L-1, 1); data_time];
end

x = [x; zeros(L-1, 1)];

x = upsample(x, lambda);

f = (fs*(-length(x)/2:length(x)/2 - 1)/length(x));
figure(1)
plot(f, abs(x));

filter = rcosine(1, lambda, 'sqrt', beta, delay);
x_filtered = conv(x, filter);


f_filtered = (fs*(-length(x_filtered)/2:length(x_filtered)/2 - 1)/length(x_filtered));
figure(2)
plot(f_filtered, abs(x_filtered));


x_pass_band = zeros(length(x_filtered), 1);
for i = 1:length(x_filtered)
    x_pass_band(i) = real(x_filtered(i)*exp(j*2*pi*fc*i*ts));
end

f0 = fc - 4e3;
f1 = fc + 4e3;
t1 = .05;
b = (f1-f0)/t1;
t = linspace(0, t1, b);
f_t = f0+b.*t;
chirp = cos(2*pi.*f_t.*t).';
z = zeros(640e3, 1);
LFM = [];
for i = 1:4
   LFM = [LFM; chirp; z]; 
end

lfm_symbols = [LFM; x_pass_band; LFM];



% Validation
lfm_symbols_base_band_I = zeros(length(lfm_symbols), 1);
lfm_symbols_base_band_Q = zeros(length(lfm_symbols), 1);
for n = 1: length(lfm_symbols)
    lfm_symbols_base_band_I(n) = lfm_symbols(n)*cos(2*pi*fc*n*ts); 
    lfm_symbols_base_band_Q(n) = -1*lfm_symbols(n)*sin(2*pi*fc*n*ts); 
end

lfm_symbols_base_band = lfm_symbols_base_band_I + j*lfm_symbols_base_band_Q;


x_ovf = conv(lfm_symbols_base_band, filter);

x_ovf_down_sampled = downsample(x_ovf(2*delay:length(x_ovf)), lambda); 

result = fft(x_ovf_down_sampled);

