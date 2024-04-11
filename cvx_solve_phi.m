function Phi = cvx_solve_phi(Sigma, U, RM, Phi)
% 使用CVX求解RIS的相位移动矩阵Θ
% 输入:
%   Sigma: 根据F和uk计算得到的矩阵Σ
%   U: 根据F和uk计算得到的矩阵U
%   RM: RIS元素的总数
% 输出:
%   Phi: RIS的相位移动向量

cvx_begin quiet
    cvx_precision low
    variable Phi(1, RM) complex;  % 定义复数变量Phi
    minimize((Phi * Sigma * Phi') + 2 * real(Phi * U));  % 目标函数
    subject to
    for n = 1:RM
        abs(Phi(n))<=1;
    end
cvx_end

end
