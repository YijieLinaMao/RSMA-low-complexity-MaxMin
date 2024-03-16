function [pc, c, MMFrate] = multicast_wmmse(H, Pt, tolerance)

[Nt, K] = size(H);

P_c = Pt;

[U, ~, ~] = svd(H);
pc = U(:, 1) * sqrt(P_c);

MMFrate_last = -inf;
maxIter = 1000;

for n = 1:maxIter

    T_c1 = square_abs(H(:, 1)'*pc) + 1;
    T_c2 = square_abs(H(:, 2)'*pc) + 1;

    g_ck = zeros(K, 1);
    g_ck(1) = pc' * H(:, 1) / T_c1;
    g_ck(2) = pc' * H(:, 2) / T_c2;

    u_ck = zeros(K, 1);
    u_ck(1) = T_c1 / log(2);
    u_ck(2) = T_c2 / log(2);

    [pc, c, MMFrate] = multicast_wmmse_update(H, Pt, g_ck, u_ck);

    fprintf("Multicast | %3d | obj = %f | |obj - obj_last| = %f\n", n, MMFrate, abs(MMFrate - MMFrate_last));

    if abs(MMFrate - MMFrate_last) <= tolerance
        break;
    end
    MMFrate_last = MMFrate;

end
