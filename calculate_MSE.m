function MSE = calculate_MSE(H, u, F, sigma_squared, L, K, Nr, Nt)
% 计算所有用户的均方误差
% 输入:
%   H: 所有AP到所有用户的组合信道矩阵，维度为[Nr, Nt, L, K]
%   u: 所有用户的组合向量，维度为[Nr, K]
%   F: 所有AP的预编码矩阵，维度为[Nt, K, L]
%   sigma_squared: 所有用户接收到的噪声功率，维度为[K, 1]
%   L: AP总数
%   K: 用户总数
%   Nr: 用户的接收天线数
%   Nt: AP的发射天线数
% 输出:
%   MSE: 所有用户的均方误差，维度为[K, 1]

% 初始化 MSE 向量
MSE = zeros(K, 1);

% 计算每个用户的 MSE
for k = 1:K
    % 合并所有L个AP的信道矩阵
    Hk = zeros(Nr, Nt * L);
    for l = 1:L
        Hk(:, (l-1)*Nt + 1:l*Nt) = H(:, :, l, k);
    end
    
    % 合并所有L个AP的预编码向量
    Fk = zeros(Nt * L, 1);
    for l = 1:L
        Fk((l-1)*Nt + 1:l*Nt) = F(:, k, l);
    end
    
    % 计算信号项和干扰项
    signal = abs(u(:, k)' * Hk * Fk)^2;
    interference = 0;
    for i = 1:K
        if i ~= k
            Fi = zeros(Nt * L, 1);
            for l = 1:L
                Fi((l-1)*Nt + 1:l*Nt) = F(:, i, l);
            end
            interference = interference + abs(u(:, k)' * Hk * Fi)^2;
        end
    end
    
    % 计算 MSE
    MSE(k) = interference + sigma_squared(k) * norm(u(:, k))^2 - 2 * real(u(:, k)' * Hk * Fk) + 1;
end

end
