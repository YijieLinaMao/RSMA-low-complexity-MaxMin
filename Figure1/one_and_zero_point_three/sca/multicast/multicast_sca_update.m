function [c, t, pc] = ...
          multicast_sca_update(H, Pt, pcPrevious)
Nt = size(H, 1);
cvx_begin quiet
    variable t
    variable alphac1
    variable alphac2
    variable rhoc1
    variable rhoc2
    variable pc(Nt) complex
    variable c(2)
    
    maximize(t)
    subject to
        t <= c(1)
        t <= c(2)

        c(1) + c(2) <= alphac1
        c(1) + c(2) <= alphac2

        alphac1 <= log(1 + rhoc1) / log(2)
        alphac2 <= log(1 + rhoc2) / log(2)

        rhoc1 <= 2*real(pcPrevious'*H(:, 1)*H(:, 1)'*pc) - ...
            abs(H(:, 1)'*pcPrevious)^2
        rhoc2 <= 2*real(pcPrevious'*H(:, 2)*H(:, 2)'*pc) - ...
            abs(H(:, 2)'*pcPrevious)^2

        pc'*pc <= Pt

        c(1) >= 0
        c(2) >= 0
cvx_end
