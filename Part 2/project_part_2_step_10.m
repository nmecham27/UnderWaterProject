clear all;
%clf;

load("EndOfPart9.mat");
load("ofdm_map.mat");
sampling_rate_192 = 192*10^3;
ts_192 = 1/sampling_rate_192;

%Find the minimum value from the p_null matrix
[min_p_null, min_p_null_index] = min(p_null, [], "all", "linear");

%Since the index originally started at 2200 (see line 178). To get the
%index into our actual data we take the found index plus the start then 
min_p_null_index = 2200 + min_p_null_index - 1;

previous_starting_index = min_p_null_index - (K+L)*lambda;

%Loop through each ZP-OFDM symbol
for w = 1:21
    for starting_index = (previous_starting_index + (K+L)*lambda) + (-2*lambda:1:2*lambda)
        yBB_w_cfo_comp = zeros(1,(K+L)*lambda);
        for n_index = 0:1:(K+L)*lambda-1
            %CFO compensation
            yBB_w_cfo_comp(n_index+1) = yBB_filtered(n_index+starting_index+1)*exp(-1i*2*pi*eps*(n_index+starting_index+1)*ts_192);
        end

        %Down sampling
        yBB_w_down_sampled = zeros(1, K+L);
        for i = 0:1:K+L-1
           yBB_w_down_sampled(i+1) = yBB_w_cfo_comp((i*lambda)+1);
        end

        %Obtaining frequency data
        %This is the biggest problem loop for why everything takes so long
        z_w_freq_data = zeros(1,K);
        for dnull = 1:1:K
            if(ofdm_map(dnull) == 0) %found a null subcarrier index so process that data
                for j = 0:1:K+L-1
                    z_w_freq_data(dnull) = yBB_w_down_sampled(j+1)*exp(-1i*((2*pi*(dnull-1)*j)/(K)));
                end
            end
        end

        % Need to add in some check to find the starting index for this run
        % and set the previous_starting_index value equal to this
    end
end

%This only calculates the Zw for the null subcarriers right now. Probably
%best to make sure all parameters are generated right and then fix the
%frequency data generation section to get all Zw. It is just going to add a
%lot to the script.