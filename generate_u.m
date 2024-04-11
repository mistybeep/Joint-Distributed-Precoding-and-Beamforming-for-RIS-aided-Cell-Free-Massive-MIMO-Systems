function U = generate_u(H, F, sigma_squared, L, K, Nr)
% 计算所有用户的组合向量
% 输入:
%   H: 所有AP到所有用户的组合信道矩阵，维度为[Nr, Nt, L, K]
%   F: 更新后的预编码矩阵，维度为[Nt, K, L]
%   sigma_squared: 所有用户接收到的噪声功率，维度为[K, 1]
%   L: AP总数
%   K: 用户总数
%   Nr: 用户的接收天线数
% 输出:
%   U: 所有用户的组合向量，维度为[Nr, K]

% 初始化组合向量矩阵
[~,Nt,~,~] = size(H);
U = zeros(Nr, K);

% 计算每个用户的组合向量
for k = 1:K
    % 计算用户k的组合信道矩阵Hk
    Hk = zeros(Nr, L * Nt);
    for l = 1:L
        Hk(:, (l - 1) * Nt + 1:l * Nt) = H(:, :, l, k);
    end

    % 初始化Wk为0矩阵
    Wk = zeros(Nr, Nr);
    for i = 1:K
        % 初始化fi为用户i的预编码向量
        fi = zeros(L * Nt, 1);
        for l = 1:L
            fi((l - 1) * Nt + 1:l * Nt) = F(:, i, l);
        end

        % 累加用户i的信号功率
        Wk = Wk + (Hk * fi) * (Hk * fi)';
    end


    % 初始化fk为用户k的预编码向量
    fk = zeros(L * Nt, 1);
    for l = 1:L
        fk((l - 1) * Nt + 1:l * Nt) = F(:, k, l);
    end


    % 计算ak
    ak = Hk * fk;

    % 计算组合向量uk
    uk = (Wk + sigma_squared(k) * eye(Nr)) \ ak;

    % 将组合向量uk存储到矩阵U中
    U(:, k) = uk;
end
end
