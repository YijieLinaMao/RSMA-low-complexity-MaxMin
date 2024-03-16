% This function calculate the max-min rate and power allocatoin
% H is channel matrix
% P is power constraint

function [t, rate] = lowComplexity(H, P)

if H(:, 1)'*H(:, 1) < H(:, 2)'*H(:, 2)
    H = H(:, [2, 1]);
end

H_bar = zeros(size(H));
H_bar(:, 1) = H(:, 1) / norm(H(:, 1));
H_bar(:, 2) = H(:, 2) / norm(H(:, 2));

rho = 1 - abs(H_bar(:, 1)' * H_bar(:, 2)) ^ 2;

coefficient = 1 / sqrt(2*(1 + abs(H_bar(:, 1)'*H_bar(:, 2))));
pc_bar = coefficient * (H_bar(:, 1) + H_bar(:, 2)*exp(-1i*angle(H_bar(:, 1)'*H_bar(:, 2))));

rho1 = H(:, 1)' * H(:, 1) * rho;
rho2 = H(:, 2)' * H(:, 2) * rho;

rhoc1 = abs(H(:, 1)'*pc_bar) ^ 2;
rhoc2 = abs(H(:, 2)'*pc_bar) ^ 2;

Gamma = max(1/rho2 - 1/rho1, 0);

tList = [0, 1];

t = Gamma / P;
if t < 1 && t > 0
    % t + 0.001 and t - 0.001 are for jump discontinuity
    tList = [tList, t, t + 0.001, t - 0.001];
end

t = (0.5*rho2*Gamma - rhoc2*P - 1)/(rho2*P - 2*rhoc2*P) - 1/rho1/P - Gamma/P/2;
if t < 1 && t > 0
    tList = [tList, t];
end

t = rhoc2 / (rho1+rhoc2);
if t < 1 && t > 0
    tList = [tList, t];
end

t = (2*rhoc2*P - rho1*Gamma - rho2*Gamma) / (rho1 - rho2 + 2*rhoc2) / P;
if t < 1 && t > 0
    tList = [tList, t];
end

rateList = [];
for t = tList
    rate = maxMinRate(H, t, P);
    rateList = [rateList, rate];
end

[rate, optimalIndex] = max(rateList);
t = tList(optimalIndex);
