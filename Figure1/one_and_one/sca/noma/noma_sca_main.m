clear all;
clc;

load channel.mat;

snrdB_list = 5:5:30;

noma_sca_time_list = [];
noma_sca_mmf_rate_list = [];

tolerance = 1e-3;

for snrdB = snrdB_list
   
    Pt = 10 ^ (snrdB/10);
    noma_mmf_rate = 0;

    tic
    for n = 1:channel_num

        H = channel(:, :, n);
        [~, ~, mmf_rate1] = noma_sca(H, Pt, tolerance);
        [~, ~, mmf_rate2] = noma_sca(H(:, [2, 1]), Pt, tolerance);
        noma_mmf_rate = noma_mmf_rate + max(mmf_rate1, mmf_rate2);

    end
    time = toc;

    noma_sca_time_list = [noma_sca_time_list, time / channel_num];
    noma_sca_mmf_rate_list = [noma_sca_mmf_rate_list, noma_mmf_rate / channel_num];

end

filename = 'noma_sca.mat';
save(filename, ...
    'K', ...
    'Nt', ...
    'channel_num', ...
    'Sigma', ...
    'snrdB_list', ...
    'noma_sca_time_list', ...
    'noma_sca_mmf_rate_list');
