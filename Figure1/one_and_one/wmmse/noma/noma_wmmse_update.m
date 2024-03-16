function [p1, p2, MMFrate] = noma_wmmse_update(H, Pt, g_11, g_21, g_22, u_11, u_21, u_22)

[Nt, K] = size(H);

cvx_begin quiet

    variable t
    variable p1(Nt) complex
    variable p2(Nt) complex
    
    T_11 = square_abs(H(:, 1)'*p1) + square_abs(H(:, 1)'*p2) + 1;
    T_21 = square_abs(H(:, 2)'*p1) + square_abs(H(:, 2)'*p2) + 1;
    T_22 = square_abs(H(:, 2)'*p2) + 1;

    epsilon_11 = square_abs(g_11)*T_11 - 2*real(g_11*H(:, 1)'*p1) + 1;
    epsilon_21 = square_abs(g_21)*T_21 - 2*real(g_21*H(:, 2)'*p1) + 1;
    epsilon_22 = square_abs(g_22)*T_22 - 2*real(g_22*H(:, 2)'*p2) + 1;

    xi_11 = u_11*epsilon_11 - log2(u_11);
    xi_21 = u_21*epsilon_21 - log2(u_21);
    xi_22 = u_22*epsilon_22 - log2(u_22);

    minimise(t)
    subject to
        t >= xi_11
        t >= xi_21
        t >= xi_22
        p1'*p1 + p2'*p2 <= Pt

cvx_end

MMFrate = 1/log(2) + log2(log(2)) - t;
