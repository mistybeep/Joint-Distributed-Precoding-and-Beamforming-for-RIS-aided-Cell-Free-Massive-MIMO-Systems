function [F, lambda] = generate_f(H, u, F, omega, Pl_max, L, K, Nt)
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
lambda = zeros(L, 1);

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

% 使用二分法求解每个AP的对偶变量lambda
for l = 1:L
    lambda_min = 0;
    lambda_max = 1;
    max_iter = 100;
    iter = 0;
    while true
        lambda(l) = (lambda_max + lambda_min) / 2;
        P_temp = 0;
        for k = 1:K
            % 提取针对AP l 的 Lambda 分块
            Lambda_ll = Lambda((l-1)*Nt+1:l*Nt, (l-1)*Nt+1:l*Nt);
            
            % 计算sum_term
            sum_term = zeros(Nt, 1);
            for ell = setdiff(1:L, l)
                Lambda_le = Lambda((l-1)*Nt+1:l*Nt, (ell-1)*Nt+1:ell*Nt);
                f_ellk = F(:, k, ell);
                sum_term = sum_term + Lambda_le * f_ellk;
            end
            
            % 计算AP l 对用户 k 的预编码向量 F_temp
            blk = b((l-1)*Nt+1:l*Nt, k);
            F_temp = (Lambda_ll + lambda(l) * eye(Nt)) \ (omega(k) * blk - sum_term);
            
            % 累加计算当前lambda下的发射功率
            P_temp = P_temp + trace(F_temp*F_temp');
        end
        
        % 调整lambda的搜索区间
        if P_temp > Pl_max
            lambda_min = lambda(l);
        else
            lambda_max = lambda(l);
        end
        
        % 迭代次数加一
        iter = iter + 1;
        
        % 检查是否达到收敛条件或最大迭代次数
        if abs(lambda_max - lambda_min) < 1e-5 || iter > max_iter
            break
        end
    end
    for k = 1:K
        % 提取Lambda的分块
        Lambda_ll = Lambda((l-1)*Nt+1:l*Nt, (l-1)*Nt+1:l*Nt);
        
        % 提取b向量的相应分量
        blk = b((l-1)*Nt+1:l*Nt, k);
        
        % 计算sum_term
        sum_term = zeros(Nt, 1);
        for ell = setdiff(1:L, l)
            Lambda_le = Lambda((l-1)*Nt+1:l*Nt, (ell-1)*Nt+1:ell*Nt);
            f_ellk = F(:, k, ell);
            sum_term = sum_term + Lambda_le * f_ellk;
        end
        
        % 更新F的相应分量
        F(:, k, l) = (Lambda_ll + lambda(l) * eye(Nt)) \ (omega(k) * blk - sum_term);
    end
end



% % 更新预编码向量
% for l = 1:L
%     for k = 1:K
%         % 提取Lambda的分块
%         Lambda_ll = Lambda((l-1)*Nt+1:l*Nt, (l-1)*Nt+1:l*Nt);
% 
%         % 提取b向量的相应分量
%         blk = b((l-1)*Nt+1:l*Nt, k);
% 
%         % 计算sum_term
%         sum_term = zeros(Nt, 1);
%         for ell = setdiff(1:L, l)
%             Lambda_le = Lambda((l-1)*Nt+1:l*Nt, (ell-1)*Nt+1:ell*Nt);
%             f_ellk = F(:, k, ell);
%             sum_term = sum_term + Lambda_le * f_ellk;
%         end
% 
%         % 更新F的相应分量
%         F(:, k, l) = (Lambda_ll + lambda(l) * eye(Nt)) \ (omega(k) * blk - sum_term);
%     end
% end

end
