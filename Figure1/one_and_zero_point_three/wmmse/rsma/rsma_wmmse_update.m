function [p1, p2, pc, c, MMFrate] = rsma_wmmse_update(H, Pt, g_ck, g_k, u_ck, u_k)

[Nt, ~] = size(H);
cvx_begin quiet
    variable t
    variable x(2)
    variable p1(Nt) complex
    variable p2(Nt) complex
    variable pc(Nt) complex

    T_c1 = square_abs(H(:, 1)'*pc) + square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_c2 = square_abs(H(:, 2)'*pc) + square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;

    epsilon_c1 = square_abs(g_ck(1))*T_c1 - 2*real(g_ck(1)*H(:, 1)'*pc) + 1;
    epsilon_c2 = square_abs(g_ck(2))*T_c2 - 2*real(g_ck(2)*H(:, 2)'*pc) + 1;

    xi_c1 = u_ck(1)*epsilon_c1 - log2(u_ck(1));
    xi_c2 = u_ck(2)*epsilon_c2 - log2(u_ck(2));

    T_1 = square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_2 = square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;

    epsilon_1 = square_abs(g_k(1))*T_1 - 2*real(g_k(1)*H(:, 1)'*p1) + 1;
    epsilon_2 = square_abs(g_k(2))*T_2 - 2*real(g_k(2)*H(:, 2)'*p2) + 1;

    xi_1 = u_k(1)*epsilon_1 - log2(u_k(1));
    xi_2 = u_k(2)*epsilon_2 - log2(u_k(2));

    minimise(t)
    subject to
        t >= xi_1 + x(1)
        t >= xi_2 + x(2)
        1/log(2) + log2(log(2)) + x(1) + x(2) >= xi_c1
        1/log(2) + log2(log(2)) + x(1) + x(2) >= xi_c2
        p1'*p1 + p2'*p2 + pc'*pc <= Pt
        x <= 0

cvx_end

MMFrate = 1/log(2) + log2(log(2)) - t;
c = -x;