function [Sigma, U] = generate_Sigma_U(Hd, Hr, G, u, F, omega, L, K, R, M)
% 生成矩阵Sigma和U
% 输入:
%   Hd: 直接信道矩阵，维度为[Nr, Nt, L, K]
%   Hr: RIS到用户的反射信道矩阵，维度为[M, Nr, R, K]
%   G: AP到RIS的信道矩阵，维度为[M, Nt, R, L]
%   u: 所有用户的组合向量，维度为[Nr, K]
%   F: 所有AP的预编码矩阵，维度为[Nt, K, L]
%   omega: 用户权重，维度为[K, 1]
%   L: AP总数
%   K: 用户总数
%   R: RIS总数
%   M: 每个RIS的元素数
% 输出:
%   Sigma: 矩阵Sigma，维度为[RM, RM]
%   U: 矩阵U，维度为[RM, 1]

Nt = size(Hd, 2);  % AP的发射天线数
RM = R * M;  % RIS元素的总数

% 初始化Sigma和U
Sigma = zeros(RM, RM);
U = zeros(RM, 1);

% 计算Sigma和U
for k = 1:K
    % 初始化ck为RIS反射信道矩阵与用户k组合向量的乘积
    ck = zeros(M * R, 1);
    for r = 1:R
        ck((r-1)*M + 1:r*M) = Hr(:, :, r, k) * u(:, k);
    end

    for i = 1:K
        term1 = zeros(M * R, 1);
        term2 = 0;
        for l = 1:L
            cdlk = Hd(:, :, l, k)' * u(:, k);
            % 构建gl
            gl = zeros(M * R, Nt);
            for r = 1:R
                gl((r-1)*M + 1:r*M, :) = G(:, :, r, l);
            end
            d_lk = diag(ck) * gl;
            F_li = F(:, i, l);
            Sigma = Sigma + omega(k) * (d_lk * F_li) * (d_lk * F_li)';
            term1 = term1 + d_lk * F_li;
            term2 = term2 + cdlk' * F_li;
        end
        U = U + omega(k) * term1 * term2;
        if i == k
            U = U - omega(k) * term1;
        end
    end
end

end
