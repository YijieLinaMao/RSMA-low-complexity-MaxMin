clear all;
clc;

load channel.mat;

snrdB_list = 5:5:30;

low_complexity_time_list = [];
low_complexity_mmf_rate_list = [];

for snrdB = snrdB_list
   
    Pt = 10 ^ (snrdB/10);
    low_complexity_mmf_rate = 0;

    tic
    for n = 1:channel_num

        H = channel(:, :, n);
        [~, mmf_rate] = lowComplexity(H, Pt);
        low_complexity_mmf_rate = low_complexity_mmf_rate + mmf_rate;

    end
    time = toc;

    low_complexity_time_list = [low_complexity_time_list, time / channel_num];
    low_complexity_mmf_rate_list = [low_complexity_mmf_rate_list, low_complexity_mmf_rate / channel_num];

end

filename = 'low_complexity.mat';
save(filename, ...
    'K', ...
    'Nt', ...
    'channel_num', ...
    'Sigma', ...
    'snrdB_list', ...
    'low_complexity_time_list', ...
    'low_complexity_mmf_rate_list');
