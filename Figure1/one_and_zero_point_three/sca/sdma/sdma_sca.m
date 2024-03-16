function [p1, p2, MMFrate] = sdma_sca(H, Pt, tolerance)

P = pinv(H');
p1Previous = sqrt(0.5 * Pt) * P(:, 1) / norm(P(:, 1));
p2Previous = sqrt(0.5 * Pt) * P(:, 2) / norm(P(:, 2));


betaPrevious = 1 + [abs(H(:, 1)' * p2Previous) ^ 2, ...
                    abs(H(:, 2)' * p1Previous) ^ 2];

tPrevious = -inf;

maxIter = 1000;
for iter = 1:maxIter
    [t, betaa, ...
     p1, p2] = ...
    sdma_sca_update(H, Pt, betaPrevious, ...
               p1Previous, p2Previous);

    fprintf("RSMA | %3d | obj = %f | |obj - obj_last| = %f\n", iter, t, abs(t - tPrevious));

    if(abs(t - tPrevious) <= tolerance)
        break;
    end

    tPrevious = t;
    betaPrevious = betaa;
    p1Previous = p1;
    p2Previous = p2;
end
MMFrate = t;
