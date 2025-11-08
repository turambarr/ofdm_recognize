function [Ng_hat, diagG] = paper_estimate_Ng(y, N_hat, Fs_hat, Np, Td_ns, Verbose, TopK)
% 依据论文(式(24),(10),(11))在 Fs_hat 采样域估计 Ng
% 输入:
%   y      : 复基带列向量（采样率应为 Fs_hat）
%   N_hat  : 已估的子载波数 N
%   Fs_hat : 已估的信道带宽 Fs（采样率）
%   Np     : 谐波数上限（建议 1~4，论文指出 Np > N/Ng 无收益）
%   Td_ns  : Ku波段延迟扩展量级，默认 108 ns
%   Verbose: 显示候选评分摘要
%   TopK   : Verbose 时显示前 K 个得分
%
% 输出:
%   Ng_hat : 循环前缀样本数估计
%   diagG  : 诊断结构（候选 xi 与打分）

    if nargin<4 || isempty(Np),    Np = 4;    end
    if nargin<5 || isempty(Td_ns), Td_ns = 108; end
    if nargin<6 || isempty(Verbose), Verbose = false; end
    if nargin<7 || isempty(TopK),    TopK = 10; end

    y = y(:);
    L = numel(y);
    if L <= N_hat + 10
        error('paper_estimate_Ng: 数据过短(<=N_hat)，无法估计 Ng');
    end

    % ---- 预计算 R^0_y(N) 所需乘积序列（式(11)的内核）----
    M = L - N_hat;                              % 可用样本数
    c = y(1+N_hat : N_hat+M) .* conj(y(1:M));   % c(n)=y(n+N)*conj(y(n))
    n = (0:M-1).';                               % 索引 n

    % ---- 候选集：xi = N + Ng，Ng 基于 Td_ns 的范围（式(24)）----
    Ng_lo = max(2, ceil( (Td_ns*1e-9 * Fs_hat) / 2 ));   % ~ Td/2 * Fs
    Ng_hi =          floor( (2*Td_ns*1e-9 * Fs_hat)    );% ~ 2Td  * Fs
    xi_list = (N_hat + Ng_lo : N_hat + Ng_hi).';         % xi = N+Ng

    % ---- 对每个 xi 累加 |R^α_y(N)| 在 |p|<=Np 的和，取最大（式(10)）----
    score = zeros(numel(xi_list),1);
    for t = 1:numel(xi_list)
        xi = xi_list(t);
        s = 0;
        for p = -Np:Np
            alpha = p / xi;
            w = exp(-1j*2*pi*alpha*n);
            Ralpha = (1/M) * (w.' * c);
            s = s + abs(Ralpha);
        end
        score(t) = s;
    end

    % 选最大得分处的 xi ⇒ Ng_hat = xi - N_hat
    [~,ix] = max(score);
    xi_star = xi_list(ix);
    Ng_hat  = xi_star - N_hat;

    if Verbose
        fprintf('[Ng] xi candidates: %d (Ng in [%d,%d], N=%d)\n', numel(xi_list), Ng_lo, Ng_hi, N_hat);
        [score_sorted, order] = sort(score, 'descend');
        K = min(TopK, numel(order));
        for k=1:K
            xi_k = xi_list(order(k)); Ng_k = xi_k - N_hat;
            mark = ''; if xi_k==xi_star, mark=' <= best'; end
            fprintf('  [Ng] #%02d: Ng=%d (xi=%d) score=%.6g%s\n', k, Ng_k, xi_k, score_sorted(k), mark);
        end
    end

    diagG = struct('xi_list',xi_list,'score',score,'Ng_lo',Ng_lo,'Ng_hi',Ng_hi,'M',M,'Np',Np);
end
