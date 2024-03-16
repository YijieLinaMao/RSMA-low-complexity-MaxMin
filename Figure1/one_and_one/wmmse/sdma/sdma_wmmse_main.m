clear all;
clc;

load channel.mat;

snrdB_list = 5:5:30;

sdma_wmmse_time_list = [];
sdma_wmmse_mmf_rate_list = [];

tolerance = 1e-3;

for snrdB = snrdB_list

    Pt = 10 ^ (snrdB/10);
    sdma_mmf_rate = 0;

    tic
    for n = 1:channel_num

        H = channel(:, :, n);
        [~, ~, mmf_rate] = sdma_wmmse(H, Pt, tolerance);
        sdma_mmf_rate = sdma_mmf_rate + mmf_rate;

    end
    time = toc;

    sdma_wmmse_time_list = [sdma_wmmse_time_list, time / channel_num];
    sdma_wmmse_mmf_rate_list = [sdma_wmmse_mmf_rate_list, sdma_mmf_rate / channel_num];

end

filename = 'sdma_wmmse.mat';
save(filename, ...
    'K', ...
    'Nt', ...
    'channel_num', ...
    'Sigma', ...
    'snrdB_list', ...
    'sdma_wmmse_time_list', ...
    'sdma_wmmse_mmf_rate_list');
