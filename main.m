clear
clc
% 初始化参数
L = 5;  % AP总数
K = 4;  % 用户总数
Nt = 4;  % AP的发射天线数
Nr = 2;  % 用户的接收天线数
R = 2;  % RIS总数
M = 100; % 每个RIS的元素数
Pl_max = 0; % AP的最大发射功率 (dBm)
Pl_max = db2pow(Pl_max);
sigma_k_squared = 10^(-80/10)*ones(K,1); % 用户接收到的噪声功率
omega = ones(K,1);
epsilon = 1e-3; % 收敛条件

% 位置设置
pos_AP = [40*(0:L-1); -50*ones(1,L); 3*ones(1,L)]';  % 所有AP的位置
pos_RIS = [60, 10, 6; 100, 10, 6];  % 两个RIS的位置
pos_user = [30, 0, 0];  % 用户的位置

% 大尺度衰落参数
C0 = -30;    % 参考距离处的路径损耗
C0 = db2pow(C0);
d0 = 1;    % 参考距离
kappa = 3; % 路径损耗指数

% 小尺度衰落参数
omega_Bu = 0.01;  % Rice因子

% 生成信道矩阵
[Hd, Hr, G] = generate_channel(Nr, Nt, L, K, R, M, pos_AP, pos_RIS, pos_user, C0, d0, kappa, omega_Bu);

% 初始化变量
F = randn(Nt, K, L) + 1i * randn(Nt, K, L);  % 所有AP的预编码矩阵
% 按功率约束调整预编码矩阵
for l = 1:L
    % 计算AP l的预编码矩阵的总功率
    Pl = 0;
    for k = 1:K
        Pl = Pl + norm(F(:, k, l), 'fro')^2;
    end

    % 如果总功率超过Pl_max，则按比例缩放预编码矩阵
    if Pl > Pl_max
        scaling_factor = sqrt(Pl_max / Pl);
        F(:, :, l) = F(:, :, l) * scaling_factor;
    end
end
trace(F(:, :, 1)*F(:, :, 1)')
Phi = exp(1i * 2 * pi * rand(1, R * M)); % 初始化RIS的相位移动向量

H = zeros(Nr, Nt, L, K);
for k = 1:K
    for l = 1:L
        H_combined = zeros(Nr, Nt);
        for r = 1:R
            % 计算每个RIS的反射信道对用户k的等效信道矩阵
            Hr_eff = Hr(:, :, r, k)' *diag(Phi((r-1)*M+1:r*M)) * G(:, :, r, l);
            % 累加所有RIS的贡献
            H_combined = H_combined +  Hr_eff;
        end
        % 将直接信道和所有RIS的反射信道相加得到组合信道
        H(:, :, l, k) = Hd(:, :, l, k) + H_combined;
    end
end

% 求初始化发射波束V后求系统和速率
u = generate_u(H, F, sigma_k_squared, L, K, Nr);
MSE = calculate_MSE(H, u, F, sigma_k_squared, L, K, Nr, Nt);
rate2 = []; % 初始化一个空向量记录rate
Rate_2 = calculate_rate(H, u, F, sigma_k_squared, K, Nr);
% Rate_2 = log2(1./MSE);
rate_old = sum(omega.*Rate_2);
rate2 = [rate2 rate_old];

% 迭代过程
max_iter = 30; % 最大迭代次数
iter1 = 1;
while(1)
    % 更新所有AP的预编码矩阵
    [F, lambda] = generate_f(H, u, F, omega, Pl_max, L, K, Nt);
    for l = 1:L
        % 计算AP l的预编码矩阵的总功率
        Pl = 0;
        for k = 1:K
            Pl = Pl + norm(F(:, k, l), 'fro')^2;
        end

        % 如果总功率超过Pl_max，则按比例缩放预编码矩阵
        if Pl > Pl_max
            scaling_factor = sqrt(Pl_max / Pl);
            F(:, :, l) = F(:, :, l) * scaling_factor;
        end
    end

    % 更新所有用户的组合向量
    u = generate_u(H, F, sigma_k_squared, L, K, Nr);

    % 生成矩阵Sigma和U
    [Sigma, U] = generate_Sigma_U(Hd, Hr, G, u, F, omega, L, K, R, M);

    % 使用CVX求解RIS的相位移动向量
    Phi_new = cvx_solve_phi(Sigma, U, R*M, Phi);
    Phi = Phi_new; % 更新RIS的相位移动向量

    for k = 1:K
        for l = 1:L
            H_combined = zeros(Nr, Nt);
            for r = 1:R
                % 计算每个RIS的反射信道对用户k的等效信道矩阵
                Hr_eff = diag(Phi((r-1)*M+1:r*M)) * G(:, :, r, l);
                % 累加所有RIS的贡献
                H_combined = H_combined + Hr(:, :, r, k)' * Hr_eff;
            end
            % 将直接信道和所有RIS的反射信道相加得到组合信道
            H(:, :, l, k) = Hd(:, :, l, k) + H_combined;
        end
    end

    rate_new = calculate_rate(H, u, F, sigma_k_squared, K, Nr);
    MSE = calculate_MSE(H, u, F, sigma_k_squared, L, K, Nr, Nt);
    % rate_new = log2(1./MSE);
    rate_new = sum(omega.*rate_new);
    rate2 = [rate2 rate_new];

    iter1 = iter1 + 1;
    if iter1 > max_iter % abs(rate_new-rate_old) / rate_old < epsilon ||
        break;
    end
    rate_old = rate_new;

end

for l = 1:L
    trace(F(:, :, l)*F(:, :, l)')
end

plot(0:iter1-1,rate2,'r-o')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')



% rate = sum_rate(H, F, u, sigma_k_squared, omega);

