function [pc, c, MMFrate] = multicast_wmmse_update(H, Pt, g_ck, u_ck)

[Nt, ~] = size(H);
cvx_begin quiet
    variable t
    variable x(2)
    variable pc(Nt) complex

    T_c1 = square_abs(H(:, 1)'*pc) + 1;
    T_c2 = square_abs(H(:, 2)'*pc) + 1;

    epsilon_c1 = square_abs(g_ck(1))*T_c1 - 2*real(g_ck(1)*H(:, 1)'*pc) + 1;
    epsilon_c2 = square_abs(g_ck(2))*T_c2 - 2*real(g_ck(2)*H(:, 2)'*pc) + 1;

    xi_c1 = u_ck(1)*epsilon_c1 - log2(u_ck(1));
    xi_c2 = u_ck(2)*epsilon_c2 - log2(u_ck(2));

    minimise(t)
    subject to
        t >= x(1)
        t >= x(2)
        1/log(2) + log2(log(2)) + x(1) + x(2) >= xi_c1
        1/log(2) + log2(log(2)) + x(1) + x(2) >= xi_c2
        pc'*pc <= Pt
        x <= 0

cvx_end

MMFrate = - t;
c = -x;
