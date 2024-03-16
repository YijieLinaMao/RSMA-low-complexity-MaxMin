% This is a code package related to the following paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F. Luo and Y. Mao, "A Practical Max-Min Fair Resource Allocation Algorithm
% for Rate-Splitting Multiple Access," in IEEE Communications Letters,
% vol. 27, no. 12, pp. 3285-3289, Dec. 2023, doi: 10.1109/LCOMM.2023.3329149.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The code is implemented in MATLAB environment with CVX toolbox 
% assisted. 
%
% Fig. 2 of the above paper will be reproduced by running this MATLAB 
% script. By changing the variable 'SNRdB', you can reproduce Fig. 2(a-c).
%
% This script relies on function 'lowComplexity.m' and 'maxMinRate.m'.
% Don't forget to put 'lowComplexity.m' and 'maxMinRate.m' in the same folder.


clear all;
clc;

% signal to noise ratio in decibel
SNRdB = 10;
Pt = 10 ^ (SNRdB/10);

% these two parameters control the density of grid
gammaNumber = 400;
thetaNumber = 300;

gammaList = linspace(0.1, 1, gammaNumber);
thetaList = linspace(0.01, pi - 0.01, thetaNumber);

optimalTList = zeros(gammaNumber, thetaNumber);

for gammaIndex = 1:gammaNumber
    for thetaIndex = 1:thetaNumber

        % channel matrix
        H = [1 / sqrt(2), gammaList(gammaIndex) / sqrt(2)
             1 / sqrt(2), gammaList(gammaIndex) / sqrt(2) * exp(-1i*thetaList(thetaIndex))];

        % calculate the max-min rate
        [t, ~] = lowComplexity(H, Pt);

        optimalTList(gammaIndex, thetaIndex) = t;

    end
end
gammadB = 20 * log10(gammaList);
rhoList = 0.5 * (1 - cos(thetaList));
[Rho, Gamma] = meshgrid(rhoList, gammadB);

figure;
contourf(Rho, Gamma, optimalTList, 30, 'ShowText', 'off', 'LineWidth', 0.01);
xlabel('\rho', 'FontSize', 20);
ylabel('channel strength disparity \gamma_{dB}[dB]', 'FontSize', 14);
title(sprintf('(c) SNR = %ddB', SNRdB), 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'Normal');
colormap default;
c = colorbar;
c.FontSize = 14;
hold on;
set(gca, 'LooseInset', [0,0,0.02,0]);
set(gca, 'FontSize', 14);
