clear variables;
%clf;

W = 24;




load("EndOfPart9.mat");
load("ofdm_map.mat");
sampling_rate_192 = 192*10^3;
ts_192 = 1/sampling_rate_192;

n0_w_results = zeros(1, W);
eps_w_results = zeros(1, W);
z_results = zeros(K, W);

%Find the minimum value from the p_null matrix
[p_null_col_mins, p_null_col_min_row] = min(p_null);
[min_p_null_value, min_p_null_col_index] = min(p_null_col_mins);

min_p_null_row_index = p_null_col_min_row(min_p_null_col_index);

n0_starting = min_p_null_row_index + 2200;
n0_eps = (min_p_null_col_index - 1) * .1 - 2;

n0_w_results(1) = n0_starting;
eps_w_results(1) = n0_eps;

%Since the index originally started at 2200 (see line 178). To get the
%index into our actual data we take the found index plus the start then
min_p_null_index = 2200 + min_p_null_row_index - 1;

%previous_starting_index = min_p_null_index - (K+L)*lambda;




n_index = 0:1:(K+L)*lambda-1;
down_sample_index = 0:1:K+L-1;

%FFT matrix
m = [1:K].';
j = [1:K+L];
fft_matrix = exp(-1i*((2*pi*(m-1)*(j-1))/(K)));
prev_n0_w = n0_starting;

p_null_n_index = 1;
p_null_eps_index = 1;


%Loop through each ZP-OFDM symbol
for w = 2:21
    p_null = zeros(length((K+L)*lambda +(-2*lambda:1:2*lambda)), length(-2:.1:2));
    for n = prev_n0_w+(K+L)*lambda+(-2*lambda:1:2*lambda)
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
        p_null_n_index = p_null_n_index + 1;
        p_null_eps_index = 1;
    end
    
    %Find the minimum value from the p_null matrix
    [p_null_col_mins, p_null_col_min_row] = min(p_null);
    [min_p_null_value, min_p_null_col_index] = min(p_null_col_mins);
    
    min_p_null_row_index = p_null_col_min_row(min_p_null_col_index);
    
    prev_n0_w = prev_n0_w +((K+L)*lambda);
    if min_p_null_row_index <= 2*lambda
       prev_n0_w = prev_n0_w - (48 - min_p_null_row_index); 
    end
    if min_p_null_row_index > 2*lambda+1
       prev_n0_w = prev_n0_w + (min_p_null_row_index - 49); 
    end
    n0_eps = (min_p_null_col_index - 1) * .1 - 2;
    
    
    p_null_n_index = 1;
    
    n0_w_results(w) = prev_n0_w;
    eps_w_results(w) = n0_eps;
    
    
    
end

%This only calculates the Zw for the null subcarriers right now. Probably
%best to make sure all parameters are generated right and then fix the
%frequency data generation section to get all Zw. It is just going to add a
%lot to the script.