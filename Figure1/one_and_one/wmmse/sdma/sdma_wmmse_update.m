function [p1, p2, MMFrate] = sdma_wmmse_update(H, Pt, g_k, u_k)

[Nt, ~] = size(H);
cvx_begin quiet
    variable t
    variable p1(Nt) complex
    variable p2(Nt) complex

    T_1 = square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_2 = square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;

    epsilon_1 = square_abs(g_k(1))*T_1 - 2*real(g_k(1)*H(:, 1)'*p1) + 1;
    epsilon_2 = square_abs(g_k(2))*T_2 - 2*real(g_k(2)*H(:, 2)'*p2) + 1;

    xi_1 = u_k(1)*epsilon_1 - log2(u_k(1));
    xi_2 = u_k(2)*epsilon_2 - log2(u_k(2));

    minimise(t)
    subject to
        t >= xi_1
        t >= xi_2
        p1'*p1 + p2'*p2 <= Pt
        
cvx_end

MMFrate = 1/log(2) + log2(log(2)) - t;
