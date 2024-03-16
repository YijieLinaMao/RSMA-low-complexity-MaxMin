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
% This script relies on function 'maxMinRate.m'.
% Don't forget to put 'maxMinRate.m' in the same folder.


clear all;
clc;

% signal to noise ratio in decibel
SNRdB = 30;
Pt = 10 ^ (SNRdB/10);

% these two parameters control the density of grid
gammaNumber = 400;
thetaNumber = 300;

gammaList = linspace(0.1, 1, gammaNumber);
thetaList = linspace(0.01, pi - 0.01, thetaNumber);

relativeGainList = zeros(gammaNumber, thetaNumber);
relativeGainOverSDMA = zeros(gammaNumber, thetaNumber);
relativeGainOverNOMA = zeros(gammaNumber, thetaNumber);
relativeGainOverMulticast = zeros(gammaNumber, thetaNumber);

for gammaIndex = 1:gammaNumber
    for thetaIndex = 1:thetaNumber

        % channel matrix
        H = [1 / sqrt(2), gammaList(gammaIndex) / sqrt(2)
             1 / sqrt(2), gammaList(gammaIndex) / sqrt(2) * exp(-1i*thetaList(thetaIndex))];

        % normalized channel vectors
        H_bar = zeros(size(H));
        H_bar(:, 1) = H(:, 1) / norm(H(:, 1));
        H_bar(:, 2) = H(:, 2) / norm(H(:, 2));

        rho = 1 - abs(H_bar(:, 1)' * H_bar(:, 2)) ^ 2;

        % common beamforming vector's direction
        pc_bar = 1 / sqrt(2 * (1 + abs(H_bar(:, 1)' * H_bar(:, 2)))) * ...
            (H_bar(:, 1) + H_bar(:, 2) * ...
            exp(-1i * angle(H_bar(:, 1)' * H_bar(:, 2))));

        rho1 = H(:, 1)' * H(:, 1) * rho;
        rho2 = H(:, 2)' * H(:, 2) * rho;

        rhoc1 = abs(H(:, 1)' * pc_bar) ^ 2;
        rhoc2 = abs(H(:, 2)' * pc_bar) ^ 2;

        % Gamma is non-negative mathematically, the max is to correct the
        % computational error
        Gamma = max(1 / rho2 - 1 / rho1, 0);

        if Pt > Gamma
            sdmaRate = maxMinRate(H, 1, Pt);
        end
        multicastRate = maxMinRate(H, 0, Pt);

        rsmaList = [];
        nomaList = [];

        % Gamma/Pt is a non-continuous point
        t = Gamma/Pt - 0.001;
        if t < 1 && t > 0
            nomaList = [nomaList, t];
        end

        % Gamma/Pt is a non-continuous point
        t = Gamma/Pt + 0.001;
        if t < 1 && t > 0
            rsmaList = [rsmaList, t];
        end

        t = rhoc2 / (rho1+rhoc2);
        if t < 1 && t > 0 && t*Pt < Gamma
            nomaList = [nomaList, t];
        end

        t = (0.5 * rho2 * Gamma - rhoc2 * Pt - 1) / ...
            (rho2 * Pt - 2 * rhoc2 * Pt) - 1 / rho1 / Pt - Gamma / Pt / 2;
        if t < 1 && t > 0 && t*Pt > Gamma
            rsmaList = [rsmaList, t];
        end
        
        % zero is non-continuous point
        t = 0.001;
        if t*Pt < Gamma
            nomaList = [nomaList, t];
        end

        t = (2 * rhoc2 * Pt - rho1 * Gamma - rho2 * Gamma) / ...
            (rho1 - rho2 + 2 * rhoc2) / Pt;
        if t < 1 && t > 0 && t*Pt > Gamma
            rsmaList = [rsmaList, t];
        end

        nomaRateList = [];
        for t = nomaList
            nomaRate = maxMinRate(H, t, Pt);
            nomaRateList = [nomaRateList, nomaRate];
        end
        nomaRate = max(nomaRateList);

        rsmaRateList = [];
        for t = rsmaList
            rsmaRate = maxMinRate(H, t, Pt);
            rsmaRateList = [rsmaRateList, rsmaRate];
        end
        rsmaRate = max(rsmaRateList);

        if Gamma - Pt >= -0.001 * Pt
            relativeGainOverSDMA(gammaIndex, thetaIndex) = 0;
            relativeGainOverNOMA(gammaIndex, thetaIndex) = 0;
            relativeGainOverMulticast(gammaIndex, thetaIndex) = 0;
            relativeGainList(gammaIndex, thetaIndex) = 0;
        else
            if Gamma <= 0.001 * Pt
                relativeGainOverNOMA(gammaIndex, thetaIndex) = 0;
            else
                relativeGainOverNOMA(gammaIndex, thetaIndex) = max((rsmaRate-nomaRate)/nomaRate, 0);
            end
            relativeGainOverSDMA(gammaIndex, thetaIndex) = max((rsmaRate-sdmaRate)/sdmaRate, 0);
            relativeGainOverMulticast(gammaIndex, thetaIndex) = max((rsmaRate-multicastRate)/multicastRate, 0);
            num = rsmaRate - max([nomaRate, sdmaRate, multicastRate]);
            den = max([nomaRate, sdmaRate, multicastRate]);
            relativeGainList(gammaIndex, thetaIndex) = max(num/den, 0);
        end
    end
end
gammadB = 20 * log10(gammaList);
rhoList = 0.5 * (1 - cos(thetaList));
[Rho, Gamma] = meshgrid(rhoList, gammadB);

figure;
contourf(Rho, Gamma, relativeGainList, 20, 'ShowText', 'off', 'LineWidth', 0.01);
xlabel('\rho', 'FontSize', 20);
ylabel('channel strength disparity \gamma_{dB}[dB]', 'FontSize', 14);
title(sprintf('SNR = %ddB', SNRdB), 'FontSize', 14);
colormap default;
c = colorbar;
c.FontSize = 14;

% change the colorbar label to percentages
tickLabels = [];
for item = c.Ticks
    tickLabels = [tickLabels, strcat(string(100*item), '%')];
end
c.TickLabels = tickLabels;

% the percentage points in the figure
hold on;
x = rhoList([90, 150, 210]);
y = gammadB([35, 120, 200, 300]);
[X, Y] = meshgrid(x, y);
plot(X, Y, 'Marker', 'x', 'LineStyle', 'none', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'white');

% the percentages in the figure
for m = [35, 120, 200, 300]
    for n = [90, 150, 210]
        txt = sprintf('(%.0f%%, %.0f%%, %.0f%%)', 100*relativeGainOverSDMA(m, n), 100*relativeGainOverNOMA(m, n), 100*relativeGainOverMulticast(m, n));
        text(rhoList(n)-0.1, gammadB(m)+1, txt, 'Color', 'white');
    end
end
set(gca, 'LooseInset', [0,0,0.035,0]);
set(gca, 'FontSize', 14);

