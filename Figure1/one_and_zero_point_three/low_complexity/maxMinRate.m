% This function calculate the max-min rate
% H is channel matrix
% tP is the power of private streams and (1 − t)P is the power of common stream
% P is power constraint

function rate = maxMinRate(H, t, P)

% normalized channel vectors
H_bar = zeros(size(H));
H_bar(:, 1) = H(:, 1) / norm(H(:, 1));
H_bar(:, 2) = H(:, 2) / norm(H(:, 2));

rho = 1 - abs(H_bar(:, 1)'*H_bar(:, 2))^2;
        
pc_bar = 1 / sqrt(2 * (1 + abs(H_bar(:, 1)' * H_bar(:, 2)))) * ...
    (H_bar(:, 1) + H_bar(:, 2) * exp(-1i * angle(H_bar(:, 1)' * H_bar(:, 2))));
        
rho1 = H(:, 1)' * H(:, 1) * rho;
rho2 = H(:, 2)' * H(:, 2) * rho;

rhoc1 = abs(H(:, 1)'*pc_bar) ^ 2;
rhoc2 = abs(H(:, 2)'*pc_bar) ^ 2;

% Gamma is non-negative mathematically, the max is to correct the computational error
Gamma = max(1/rho2 - 1/rho1, 0);

if t*P > Gamma
    % RSMA, SDMA
    P1 = 0.5 * (t*P + Gamma);
    P2 = 0.5 * (t*P - Gamma);
    Pc = (1-t) * P;
    R1 = log2(1 + rho1*P1);
    R2 = log2(1 + rho2*P2);
    Rc = log2(1 + rhoc1*Pc/(1 + rho1*P1));
    rate = R1 + R2 + Rc - max(R1 - R2 - Rc, 0);
    rate = 0.5 * rate;
else
    P1 = t * P;
    Pc = (1-t) * P;
    if t > 1e-6
        % NOMA
        rate = log2(1 + min([rho1 * P1, rhoc2 * Pc]));
    else
        % Multicast
        Rc = log2(1 + rhoc2*Pc);
        rate = 0.5 * Rc;
    end
end