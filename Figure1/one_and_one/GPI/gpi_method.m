function [P, p_c, c, mmf_rate] = gpi_method(H, Pt, tolerance)

[Nt, K] = size(H);

noise = ones(K, 1);

P_c = Pt * 0.5;
P_p = (Pt-P_c) / K;

[U, ~, ~] = svd(H);
p_c = U(:,1) * sqrt(P_c);
P = H ./ vecnorm(H) * sqrt(P_p);

gamma_min = 0;
gamma_max = 100;
gamma_num = 1001;

alpha = 0.1;

T_k = sum(square_abs(H'*P), 2) + noise;
T_ck = T_k + square_abs(H'*p_c);

rate_c = log2(T_ck ./ T_k);
rate_p = log2(T_k ./ (T_k-square_abs(diag(H'*P))));

rate_p_sorted = sort(rate_p, 'ascend');
x = min(rate_c);
y = min((x+cumsum(rate_p_sorted)) ./ (1:K)');

max_mmf_rate = y;
corresponding_f = [p_c; reshape(P, K*Nt, 1)];

f = corresponding_f;
c = max(y - rate_p, 0);
f_last = f;

for gamma = linspace(gamma_min, gamma_max, gamma_num)
    
    for n = 1:1000

        A = zeros(Nt*(K+1), Nt*(K+1), K);
        B = zeros(Nt*(K+1));
        C = zeros(Nt*(K+1));
        D = zeros(Nt*(K+1));

        for k = 1:K
            A_k = zeros(Nt);
            for m = 1:K
                A_k = blkdiag(A_k, H(:, k)*H(:, k)');
            end
            A(:, :, k) = A_k + eye(Nt*(K+1))/Pt;
            B(:, :, k) = A(:, :, k) - blkdiag(zeros(Nt*k), H(:, k)*H(:, k)', zeros(Nt*(K-k)));
            C(:, :, k) = A(:, :, k) + blkdiag(H(:, k)*H(:, k)', zeros(Nt*K));
            D(:, :, k) = A(:, :, k);
        end

        T_k = sum(square_abs(H'*P), 2) + noise;
        T_ck = T_k + square_abs(H'*p_c);
        
        rate_c = log2(T_ck ./ T_k);
        rate_p = log2(T_k ./ (T_k-square_abs(diag(H'*P))));

        ratio_c = exp(-1/alpha*(rate_c)) / sum(exp(-1/alpha*(rate_c)));
        ratio_t = exp(-1/alpha*(c+rate_p)) / sum(exp(-1/alpha*(c+rate_p)));

        E = zeros(Nt*(K+1));
        F = zeros(Nt*(K+1));
        for k = 1:K
            E = E + ratio_t(k)*A(:, :, k)/(f'*A(:, :, k)*f) + ratio_c(k)*gamma*C(:, :, k)/(f'*C(:, :, k)*f);
            F = F + ratio_t(k)*B(:, :, k)/(f'*B(:, :, k)*f) + ratio_c(k)*gamma*D(:, :, k)/(f'*D(:, :, k)*f);
        end

        f = F \ (E*f);
        f = f / norm(f) * sqrt(Pt);

        e = norm(f - f_last, 'inf');
        if e < tolerance
            break;
        end

        f_last = f;
    end

    p_c = f(1:Nt);
    P = reshape(f(Nt+1:end), Nt, K);

    T_k = sum(square_abs(H'*P), 2) + noise;
    T_ck = T_k + square_abs(H'*p_c);
    
    rate_c = log2(T_ck ./ T_k);
    rate_p = log2(T_k ./ (T_k-square_abs(diag(H'*P))));

    rate_p_sorted = sort(rate_p, 'ascend');
    x = min(rate_c);
    y = min((x+cumsum(rate_p_sorted)) ./ (1:K)');

    c = max(y - rate_p, 0);

    if y > max_mmf_rate
        max_mmf_rate = y;
        corresponding_f = f;
    end

end

f = corresponding_f;
p_c = f(1:Nt);

P = reshape(f(Nt+1:end), Nt, K);

T_k = sum(square_abs(H'*P), 2) + noise;
T_ck = T_k + square_abs(H'*p_c);

rate_c = log2(T_ck ./ T_k);
rate_p = log2(T_k ./ (T_k-square_abs(diag(H'*P))));

rate_p_sorted = sort(rate_p, 'ascend');
x = min(rate_c);
y = min((x+cumsum(rate_p_sorted)) ./ (1:K)');
c = max(y - rate_p, 0);
mmf_rate = y;
