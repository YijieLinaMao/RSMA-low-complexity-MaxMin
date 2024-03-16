clear all;
clc;

load channel.mat;

snrdB_list = 5:5:30;

multicast_wmmse_time_list = [];
multicast_wmmse_mmf_rate_list = [];

tolerance = 1e-3;

for snrdB = snrdB_list

    Pt = 10 ^ (snrdB/10);
    multicast_mmf_rate = 0;

    tic;
    for n = 1:channel_num

        H = channel(:, :, n);
        [~, ~, mmf_rate] = multicast_wmmse(H, Pt, tolerance);
        multicast_mmf_rate = multicast_mmf_rate + mmf_rate;

    end
    time = toc;

    multicast_wmmse_time_list = [multicast_wmmse_time_list, time / channel_num];
    multicast_wmmse_mmf_rate_list = [multicast_wmmse_mmf_rate_list, multicast_mmf_rate / channel_num];

end

filename = 'multicast_wmmse.mat';
save(filename, ...
    'K', ...
    'Nt', ...
    'channel_num', ...
    'Sigma', ...
    'snrdB_list', ...
    'multicast_wmmse_time_list', ...
    'multicast_wmmse_mmf_rate_list');
