function [t, betaa, ...
          p1, p2] = ...
         sdma_sca_update(H, Pt, betaPrevious, ...
                    p1Previous, p2Previous)
Nt = size(H, 1);
cvx_begin quiet
    variable t
    variable alpha1
    variable alpha2
    variable beta1
    variable beta2
    variable rho1
    variable rho2
    variable p1(Nt) complex
    variable p2(Nt) complex

    maximize(t)
    subject to
        t <= alpha1
        t <= alpha2

        log(1 + rho1) / log(2) >= alpha1
        log(1 + rho2) / log(2) >= alpha2

        2 * real(p1Previous' * H(:, 1) * H(:, 1)' * p1) / ...
            betaPrevious(1) - ...
            square_abs(H(:, 1)' * p1Previous / ...
            betaPrevious(1)) * beta1 >= rho1
        2 * real(p2Previous' * H(:, 2) * H(:, 2)' * p2) / ...
            betaPrevious(2) - ...
            square_abs(H(:, 2)' * p2Previous / ...
            betaPrevious(2)) * beta2 >= rho2

        beta1 >= square_abs(H(:, 1)' * p2) + 1
        beta2 >= square_abs(H(:, 2)' * p1) + 1

        p1' * p1 + p2' * p2 <= Pt
cvx_end

betaa = [beta1, beta2];
