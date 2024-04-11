function R = calculate_rate(H, u, F, sigma_squared, K, Nr)
% 计算所有用户的速率
% 输入:
%   H: 所有AP到所有用户的组合信道矩阵，维度为[Nr, Nt, L, K]
%   u: 所有用户的组合向量，维度为[Nr, K]
%   F: 所有AP的预编码矩阵，维度为[Nt, K, L]
%   sigma_squared: 所有用户接收到的噪声功率，维度为[K, 1]
%   K: 用户总数
%   Nr: 用户的接收天线数
% 输出:
%   R: 所有用户的速率，维度为[K, 1]

% 初始化速率向量
R = zeros(K, 1);
[~,Nt,L,~] = size(H);
% 计算每个用户的速率
for k = 1:K
    % 计算用户k的组合信道矩阵Hk
    Hk = reshape(H(:, :, :, k), [Nr, L * Nt]);

    % 计算用户k的预编码向量fk
    fk = reshape(F(:, k, :), [L * Nt, 1]);

    % 计算信号项和干扰项
    signal = abs(u(:, k)' * Hk * fk)^2;
    interference = 0;
    for i = 1:K
        if i ~= k
            fi = reshape(F(:, i, :), [L * Nt, 1]);
            interference = interference + abs(u(:, k)' * Hk * fi)^2;
        end
    end

    % 计算SINRk
    SINRk = signal / (interference + sigma_squared(k) * norm(u(:, k))^2);

    % 计算速率Rk
    R(k) = log2(1 + SINRk);
end
end
