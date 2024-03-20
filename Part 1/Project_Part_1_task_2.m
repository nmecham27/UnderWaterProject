M = 2048; % # of subcarriers
L = 5; % # of channel taps
N = 16; % # of OFDM symbols
lambda = 24; % oversampling factor
fs_fd = 24;
beta = .125; % roll-off factor
delay = 100; %unknown delay
h = [0.227 0.46 0.688 0.46 0.227];

final_symbols = [];
starting_data = [];
total_ofdm_symbols = [];

for j = 1:N

    % Generate M QPSK symbols
    data = randi([0 1], 2*M, 1);
    data_freq = nrSymbolModulate(data, 'QPSK');
    
    starting_data = [starting_data; data_freq];

    %figure; scatter(real(data_freq), imag(data_freq));
    
    % Convert to time domaindata_freq
    data_time = ifft(data_freq);

    % Add the zero padding
    total_ofdm_symbols = [total_ofdm_symbols; zeros(L-1, 1); data_time];
end

total_ofdm_symbols = [total_ofdm_symbols; zeros(L-1,1)];
total_ofdm_symbols = total_ofdm_symbols.';


total_ofdm_symbols_upsampled = upsample(total_ofdm_symbols, lambda);

rcos_filter = rcosine(1,fs_fd,'sqrt', beta, delay);

pre_filtered_ofdm = fftshift(fft((total_ofdm_symbols_upsampled)));

filtered_ofdm_symbols = conv(total_ofdm_symbols_upsampled,rcos_filter);

post_filtered_ofdm = fftshift(fft((filtered_ofdm_symbols)));


%Trying to figure out the plotting still
figure(1);
x_plot = 1:length(pre_filtered_ofdm);
plot(x_plot, pre_filtered_ofdm);



% plot(post_filtered_ofdm, length(post_filtered_ofdm));

%     % Transmission
%     y_time = zeros(M+2*(L-1),1);
%     for m = L:M+2*(L-1)
%         for i = 1:length(h)
%             y_time(m) = y_time(m)+h(i)*x(m-i+1);
%         end
%     end
%     
%     % concat all the ofdm symbols together into an array
%     temp =  y_time(L:length(y_time)).';
%     total_ofdm_symbols = [total_ofdm_symbols, temp];
    
%     % .' does the transpose. Just ' is the hermitian
%     temp =  y_time(L:length(y_time)).';
%     
% 
%     y_time = y_time(L:length(y_time)).';
%     y_time(1:L-1) = y_time(1:L-1) + y_time(length(y_time)-(L-2):length(y_time));
%     y_time = y_time(1:length(y_time)-(L-1));
%     
%     
%     y_freq = fft(y_time);
%     h_freq = fft(h,M);
%     
%     result = y_freq./h_freq;
%     
%     final_symbols = [final_symbols; result.'];
% 
%     %figure; scatter(real(final_symbols), imag(final_symbols));





