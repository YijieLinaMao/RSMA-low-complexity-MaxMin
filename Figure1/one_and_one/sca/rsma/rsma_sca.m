function [p1, p2, pc, c, MMFrate] = rsma_sca(H, Pt, tolerance)

[Nt, K] = size(H);

P_c = Pt * 0.5;
P_p = (Pt-P_c) / K;

[U, ~, ~] = svd(H);
pcPrevious = U(:, 1) * sqrt(P_c);
P = H ./ vecnorm(H) * sqrt(P_p);
p1Previous = P(:, 1);
p2Previous = P(:, 2);

betaPrevious = 1 + [abs(H(:, 1)' * p2Previous) ^ 2, ...
                    abs(H(:, 2)' * p1Previous) ^ 2];
betacPrevious = 1 + [abs(H(:, 1)' * p1Previous) ^ 2 + ...
                     abs(H(:, 1)' * p2Previous) ^ 2, ...
                     abs(H(:, 2)' * p1Previous) ^ 2 + ...
                     abs(H(:, 2)' * p2Previous) ^ 2];

tPrevious = -inf;

maxIter = 1000;
for iter = 1:maxIter
    
    [c, t, ...
     betaa, betac, ...
     p1, p2, pc] = ...
    rsma_sca_update(H, Pt, ...
                    betaPrevious, betacPrevious, ...
                    p1Previous, p2Previous, pcPrevious);

    fprintf("RSMA | %3d | obj = %f | |obj - obj_last| = %f\n", iter, t, abs(t - tPrevious));

    if(abs(t - tPrevious) <= tolerance)
        break;
    end

    tPrevious = t;
    betaPrevious = betaa;
    betacPrevious = betac;
    p1Previous = p1;
    p2Previous = p2;
    pcPrevious = pc;
end
MMFrate = t;
