function [p1, p2, MMFrate] = sdma_wmmse(H, Pt, tolerance)

[Nt, K] = size(H);

P_p = Pt / K;

P = H ./ vecnorm(H) * sqrt(P_p);
p1 = P(:, 1);
p2 = P(:, 2);

MMFrate_last = -inf;
maxIter = 1000;

for n = 1:maxIter

    T_1 = square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_2 = square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;

    g_k = zeros(K, 1);
    g_k(1) = p1' * H(:, 1) / T_1;
    g_k(2) = p2' * H(:, 2) / T_2;

    u_k = zeros(K, 1);
    u_k(1) = T_1 / (T_1 - square_abs(H(:, 1)'*p1)) / log(2);
    u_k(2) = T_2 / (T_2 - square_abs(H(:, 2)'*p2)) / log(2);

    [p1, p2, MMFrate] = sdma_wmmse_update(H, Pt, g_k, u_k);

    fprintf("SDMA | %3d | obj = %f | |obj - obj_last| = %f\n", n, MMFrate, abs(MMFrate - MMFrate_last));

    if abs(MMFrate - MMFrate_last) <= tolerance
        break;
    end
    MMFrate_last = MMFrate;

end
