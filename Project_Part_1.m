M = 16; % # of subcarriers
L = 5; % # of channel taps
N = 4; % # of OFDM symbols
h = [0.227 0.46 0.688 0.46 0.227];

% Generate M QPSK symbols
data = randi([0 1], 2*M, 1);
data_freq = nrSymbolModulate(data, 'QPSK');

% Convert to time domain
data_time = ifft(data_freq);

% Add CP's
x = data_time(length(data_time) - L + 2: length(data_time));
x = [x; data_time];

% Transmission
y_time = zeros(M+L-1,1);
for m = L:M+L-1
    for i = 1:length(h)
        %y_time(m-L+1) = y_time(m-L+1) + h(i)*x(m-i+1);
        y_time(m) = y_time(m)+h(i)*x(m-i+1);
    end
end

y_time = y_time(L:length(y_time))';


y_freq = fft(y_time);
h_freq = fft(h,16);

result = y_freq./h_freq;







