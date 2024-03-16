function [c, t, ...
          betaa, betac, ...
          p1, p2, pc] = ...
          rsma_sca_update(H, Pt, ...
                          betaPrevious, betacPrevious, ...
                          p1Previous, p2Previous, pcPrevious)
Nt = size(H, 1);
cvx_begin quiet
    variable t
    variable alpha1
    variable alpha2
    variable alphac1
    variable alphac2
    variable beta1
    variable beta2
    variable betac1
    variable betac2
    variable rho1
    variable rho2
    variable rhoc1
    variable rhoc2
    variable p1(Nt) complex
    variable p2(Nt) complex
    variable pc(Nt) complex
    variable c(2)
    
    maximize(t)
    subject to

        c(1) >= 0
        c(2) >= 0

        c(1) + c(2) <= alphac1
        c(1) + c(2) <= alphac2

        t <= alpha1 + c(1)
        t <= alpha2 + c(2)

        log(1 + rho1) / log(2) >= alpha1
        log(1 + rho2) / log(2) >= alpha2
        log(1 + rhoc1) / log(2) >= alphac1
        log(1 + rhoc2) / log(2) >= alphac2

        2 * real(p1Previous' * H(:, 1) * H(:, 1)' * p1) / ...
            betaPrevious(1) - ...
            square_abs(H(:, 1)' * p1Previous / ...
            betaPrevious(1)) * beta1 >= rho1
        2 * real(p2Previous' * H(:, 2) * H(:, 2)' * p2) / ...
            betaPrevious(2) - ...
            square_abs(H(:, 2)' * p2Previous / ...
            betaPrevious(2)) * beta2 >= rho2

        2 * real(pcPrevious' * H(:, 1) * H(:, 1)' * pc) / ...
            betacPrevious(1) - ...
            square_abs(H(:, 1)' * pcPrevious / ...
            betacPrevious(1)) * betac1 >= rhoc1
        2 * real(pcPrevious' * H(:, 2) * H(:, 2)' * pc) / ...
            betacPrevious(2) - ...
            square_abs(H(:, 2)' * pcPrevious / ...
            betacPrevious(2)) * betac2 >= rhoc2

        beta1 >= square_abs(H(:, 1)' * p2) + 1
        beta2 >= square_abs(H(:, 2)' * p1) + 1

        betac1 >= square_abs(H(:, 1)' * p1) + square_abs(H(:, 1)' * p2) + 1
        betac2 >= square_abs(H(:, 2)' * p1) + square_abs(H(:, 2)' * p2) + 1
        
        p1' * p1 + p2' * p2 + pc' * pc <= Pt
cvx_end

betaa = [beta1, beta2];
betac = [betac1, betac2];
