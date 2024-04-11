function [Hd, Hr, G] = generate_channel(Nr, Nt, L, K, R, M, pos_AP, pos_RIS, pos_user, C0, d0, kappa, omega_Bu)
    % 大尺度衰落计算
    function L = large_scale_fading(d)
        L = C0 * (d / d0).^(-kappa);
    end

    % 小尺度衰落计算
    function H = small_scale_fading(L, rows, cols)
        H_LOS = exp(1i * 2 * pi * rand(rows, cols));  % LoS分量
        H_NLOS = (randn(rows, cols) + 1i * randn(rows, cols)) / sqrt(2);  % NLoS分量
        H = sqrt(omega_Bu / (1 + omega_Bu)) * H_LOS + sqrt(1 / (1 + omega_Bu)) * H_NLOS;
        H = sqrt(L) * H;  % 考虑大尺度衰落
    end

    % 初始化信道矩阵
    Hd = zeros(Nr, Nt, L, K);
    Hr = zeros(M, Nr, R, K);
    G = zeros(M, Nt, R, L);

    % 生成信道矩阵
    for l = 1:L
        for k = 1:K
            % 计算距离
            d_Bu = norm(pos_AP(l, :) - pos_user);  % AP到用户的距离
            L_Bu = large_scale_fading(d_Bu);  % 大尺度衰落
            Hd(:, :, l, k) = small_scale_fading(L_Bu, Nr, Nt);  % 直接信道矩阵
        end

        for r = 1:R
            % 计算距离
            d_BR = norm(pos_AP(l, :) - pos_RIS(r, :));  % AP到RIS的距离
            L_BR = large_scale_fading(d_BR);  % 大尺度衰落
            G(:, :, r, l) = small_scale_fading(L_BR, M, Nt);  % AP到RIS的信道矩阵

            for k = 1:K
                % 计算距离
                d_Ru = norm(pos_RIS(r, :) - pos_user);  % RIS到用户的距离
                L_Ru = large_scale_fading(d_Ru);  % 大尺度衰落
                Hr(:, :, r, k) = small_scale_fading(L_Ru, M, Nr);  % RIS到用户的反射信道矩阵
            end
        end
    end
end
