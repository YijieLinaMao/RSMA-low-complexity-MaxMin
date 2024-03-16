function [p1, p2, pc, c, MMFrate] = rsma_wmmse(H, Pt, tolerance)

[Nt, K] = size(H);

P_c = Pt * 0.5;
P_p = (Pt-P_c) / K;

[U, ~, ~] = svd(H);
pc = U(:, 1) * sqrt(P_c);
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

    T_c1 = square_abs(H(:, 1)'*pc) + square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_c2 = square_abs(H(:, 2)'*pc) + square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;

    g_ck = zeros(K, 1);
    g_ck(1) = pc' * H(:, 1) / T_c1;
    g_ck(2) = pc' * H(:, 2) / T_c2;

    u_ck = zeros(K, 1);
    u_ck(1) = T_c1 / T_1 / log(2);
    u_ck(2) = T_c2 / T_2 / log(2);

    [p1, p2, pc, c, MMFrate] = rsma_wmmse_update(H, Pt, g_ck, g_k, u_ck, u_k);

    fprintf("RSMA | %3d | obj = %f | |obj - obj_last| = %f\n", n, MMFrate, abs(MMFrate - MMFrate_last));

    if abs(MMFrate - MMFrate_last) <= tolerance
        break;
    end
    MMFrate_last = MMFrate;

end
