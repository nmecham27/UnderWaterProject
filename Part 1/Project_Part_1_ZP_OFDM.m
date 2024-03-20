M = 16; % # of subcarriers
L = 5; % # of channel taps
N = 4; % # of OFDM symbols
h = [0.227 0.46 0.688 0.46 0.227];

final_symbols = [];
starting_data = [];

for j = 1:N
    % Generate M QPSK symbols
    data = randi([0 1], 2*M, 1);
    data_freq = nrSymbolModulate(data, 'QPSK');
    
    starting_data = [starting_data; data_freq];

    %figure; scatter(real(data_freq), imag(data_freq));
    
    % Convert to time domaindata_freq
    data_time = ifft(data_freq);
    
    % Add ZP's
    %x = data_time(length(data_time) - L + 2: length(data_time));
    x = zeros(L-1, 1);
    x = [x; data_time; x];
    
    % Transmission
    y_time = zeros(M+2*(L-1),1);
    for m = L:M+2*(L-1)
        for i = 1:length(h)
            y_time(m) = y_time(m)+h(i)*x(m-i+1);
        end
    end
    
    % .' does the transpose. Just ' is the hermitian
    y_time = y_time(L:length(y_time)).';
    y_time(1:L-1) = y_time(1:L-1) + y_time(length(y_time)-(L-2):length(y_time));
    y_time = y_time(1:length(y_time)-(L-1));
    
    
    y_freq = fft(y_time);
    h_freq = fft(h,16);
    
    result = y_freq./h_freq;
    
    final_symbols = [final_symbols; result.'];

    %figure; scatter(real(final_symbols), imag(final_symbols));
end





