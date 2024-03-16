clear all;
clc;

load channel.mat;

snrdB_list = 5:5:30;

tolerance = 1e-3;

gpi_method_time_list = [];
gpi_method_mmf_rate_list = [];

for snrdB = snrdB_list
   
    Pt = 10 ^ (snrdB/10);
    gpi_method_mmf_rate = 0;

    tic
    for n = 1:channel_num

        H = channel(:, :, n);
        [~, ~, ~, mmf_rate] = gpi_method(H, Pt, tolerance);
        gpi_method_mmf_rate = gpi_method_mmf_rate + mmf_rate;

    end
    time = toc;

    gpi_method_time_list = [gpi_method_time_list, time / channel_num];
    gpi_method_mmf_rate_list = [gpi_method_mmf_rate_list, gpi_method_mmf_rate / channel_num];

end

filename = 'gpi_method.mat';
save(filename, ...
    'K', ...
    'Nt', ...
    'channel_num', ...
    'Sigma', ...
    'snrdB_list', ...
    'gpi_method_time_list', ...
    'gpi_method_mmf_rate_list');
