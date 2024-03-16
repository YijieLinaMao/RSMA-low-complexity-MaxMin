function [p1, p2, MMFrate] = noma_wmmse(H, Pt, tolerance)

[Nt, K] = size(H);

P_p = Pt / K;

P = H ./ vecnorm(H) * sqrt(P_p);
p1 = P(:, 1);
p2 = P(:, 2);

MMFrate_last = -inf;
maxIter = 1000;

for n = 1:maxIter

    T_11 = square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_21 = square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;
    T_22 = square_abs(H(:, 2)'*p2) + 1;

    g_11 = p1' * H(:, 1) / T_11;
    g_21 = p1' * H(:, 2) / T_21;
    g_22 = p2' * H(:, 2) / T_22;

    u_11 = T_11 / (T_11 - square_abs(H(:, 1)'*p1)) / log(2);
    u_21 = T_21 / T_22 / log(2);
    u_22 = T_22 / log(2);

    [p1, p2, MMFrate] = noma_wmmse_update(H, Pt, g_11, g_21, g_22, u_11, u_21, u_22);
    
    fprintf("NOMA | %3d | obj = %f | |obj - obj_last| = %f\n", n, MMFrate, abs(MMFrate - MMFrate_last));

    if abs(MMFrate - MMFrate_last) <= tolerance
        break;
    end
    MMFrate_last = MMFrate;

end
