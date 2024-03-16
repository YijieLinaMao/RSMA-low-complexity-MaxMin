function [p1, p2, MMFrate] = noma(H, Pt, tolerance)

[Nt, K] = size(H);

P_p = Pt / K;

P = H ./ vecnorm(H) * sqrt(P_p);
p1Previous = P(:, 1);
p2Previous = P(:, 2);

betaPrevious = 1 + [abs(H(:, 1)' * p1Previous) ^ 2, ...
                    abs(H(:, 2)' * p1Previous) ^ 2];

tPrevious = -inf;

maxIter = 1000;
for iter = 1:maxIter
    
    [t, betaa, ...
     p1, p2] = ...
    nomaUpdate(H, Pt, betaPrevious, ...
               p1Previous, p2Previous);

    fprintf("NOMA | %3d | obj = %f | |obj - obj_last| = %f\n", iter, t, abs(t - tPrevious));

    if(abs(t - tPrevious) <= tolerance)
        break;
    end

    tPrevious = t;
    betaPrevious = betaa;
    p1Previous = p1;
    p2Previous = p2;
end
MMFrate = t;
