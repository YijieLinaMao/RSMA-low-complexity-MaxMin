% This script generates 100 channel random realizations.
clear all;
clc;
K = 2;
Nt = 2;
channel_num = 100;
Sigma = [1, 1];

channel = generate_channel(Nt, K, channel_num, Sigma);
filename = 'channel.mat';

save(filename);

function channel = generate_channel(Nt, K, channel_num, Sigma)
channel = sqrt(Sigma/2) .* (randn(Nt, K, channel_num) + 1j*randn(Nt, K, channel_num));
end
