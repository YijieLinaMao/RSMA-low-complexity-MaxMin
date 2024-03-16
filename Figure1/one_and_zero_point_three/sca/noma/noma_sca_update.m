function [t, betaa, ...
          p1, p2] = ...
         noma_sca_update(H, Pt, betaPrevious, ...
                    p1Previous, p2Previous)
Nt = size(H, 1);
cvx_begin quiet
    variable t
    variable alpha1
    variable alpha2
    variable beta12
    variable beta22
    variable rho1
    variable rho12
    variable rho22
    variable p1(Nt) complex
    variable p2(Nt) complex
    
    maximize(t)
    subject to

        t <= alpha1
        t <= alpha2

        alpha1 <= log(1 + rho1) / log(2)
        alpha2 <= log(1 + rho12) / log(2)
        alpha2 <= log(1 + rho22) / log(2)

        rho1 <= 2 * real(p1Previous' * H(:, 1) * H(:, 1)' * p1) - ...
            abs(H(:, 1)' * p1Previous) ^ 2

        2 * real(p2Previous' * H(:, 1) * H(:, 1)' * p2) / ...
            betaPrevious(1) - ...
            square_abs(H(:, 1)' * p2Previous / ...
            betaPrevious(1)) * beta12 >= rho12
        2 * real(p2Previous' * H(:, 2) * H(:, 2)' * p2) / ...
            betaPrevious(2) - ...
            square_abs(H(:, 2)' * p2Previous / ...
            betaPrevious(2)) * beta22 >= rho22

        beta12 >= square_abs(H(:, 1)' * p1) + 1
        beta22 >= square_abs(H(:, 2)' * p1) + 1
        
        p1' * p1 + p2' * p2 <= Pt
cvx_end

betaa = [beta12, beta22];
