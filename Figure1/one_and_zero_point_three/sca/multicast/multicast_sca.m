function [c, pc, MMFrate] = multicast_sca(H, Pt, tolerance)

[Nt, K] = size(H);

P_c = Pt;

[U, ~, ~] = svd(H);
pcPrevious = U(:, 1) * sqrt(P_c);

tPrevious = -inf;

maxIter = 1000;
for iter = 1:maxIter

    [c, t, pc] = multicast_sca_update(H, Pt, pcPrevious);
    
    fprintf("Multicast | %3d | obj = %f | |obj - obj_last| = %f\n", iter, t, abs(t - tPrevious));

    if(abs(t - tPrevious) <= tolerance)
        break;
    end

    tPrevious = t;
    pcPrevious = pc;
end
MMFrate = t;
