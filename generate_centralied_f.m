function [F, lambda] = generate_centralied_f(H, u, F, omega, Pl_max, L, K, Nt)
% 更新预编码向量的函数，同时求解对偶变量lambda
% 输入:
%   H: 所有AP到所有用户的组合信道矩阵，维度为[Nr, Nt, L, K]
%   u: 所有用户的组合向量，维度为[Nr, K]
%   omega: 用户权重，维度为[K, 1]
%   Pl_max: AP的最大发射功率
%   L: AP总数
%   K: 用户总数
%   Nt: AP的发射天线数
% 输出:
%   F: 更新后的预编码矩阵，维度为[Nt, K, L]
%   lambda: 求解得到的对偶变量，维度为[L, 1]

% 初始化预编码矩阵和对偶变量lambda
w = zeros(Nt * L, K);
lambda = 0.5*ones(L, 1);

% 初始化 Lambda 矩阵和 b 向量
Lambda = zeros(Nt * L, Nt * L);
b = zeros(Nt * L, K);


for k = 1:K
    % 对于用户 k，首先计算串联向量 bk
    bk = [];
    for l = 1:L
        Hkl = H(:, :, l, k);  % Hkl 是 Nr x Nt 矩阵
        blk = Hkl' * u(:, k);  % blk 是 Nt x 1 向量
        bk = [bk; blk];        % bk 是 Nt*L x 1 向量
    end
    b(:, k) = bk;

    % 累加计算 Lambda 矩阵
    Lambda = Lambda + omega(k) * (bk * bk');
end

% 椭球法迭代更新lambda
epsilon = 1e-3*ones(L,1);  % 终止条件
max_iter = 10000;  % 最大迭代次数
iter = 0;
while iter < max_iter
    iter = iter + 1;
    % 更新预编码向量和计算总功率
    total_power = zeros(L,1);
    for k = 1:K
        term = zeros(Nt * L, Nt * L);
        for l=1:L
            El = [zeros(Nt,(l-1)*Nt) eye(Nt) zeros(Nt,(L-l)*Nt)];
            term = term + lambda(l)*El'*El;
        end
        w(:,k) = (Lambda + term) \ b(:, k) * omega(k);
    end
    for l=1:L
        for k=1:K
            F(:, k, l) = [zeros(Nt,(l-1)*Nt) eye(Nt) zeros(Nt,(L-l)*Nt)]*w(:,k);
            total_power(l) = total_power(l) + norm(F(:, k, l), 'fro')^2;
        end
    end
    % 检查功率约束是否满足
    if abs(Pl_max*ones(L,1)-total_power) <= epsilon
        break;  % 如果满足约束，停止迭代
    end
    % 更新lambda
    step_size = 10*epsilon .* abs(Pl_max*ones(L,1)-total_power);  % 自适应步长
    lambda = lambda + step_size .* (total_power - Pl_max*ones(L,1));
    lambda = max(lambda, 0);  % 确保lambda非负
end

end
