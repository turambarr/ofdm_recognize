function [Ng_hat, diagG] = paper_estimate_Ng(y, N_hat, Fs_hat, Np, Td_ns, Verbose, TopK, optsNg)
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
    if nargin<8 || isempty(optsNg),  optsNg = struct(); end
    % 偏好 2 的幂次（柔性选择）
    if ~isfield(optsNg,'PreferPow2'),    optsNg.PreferPow2 = false; end
    if ~isfield(optsNg,'Pow2List'),      optsNg.Pow2List   = 2.^(3:12); end  % 8,16,32,...,4096
    if ~isfield(optsNg,'SoftKeepRatio'), optsNg.SoftKeepRatio = 0.98; end    % 分数≥0.98·max 视为等优
    if ~isfield(optsNg,'MaxSnapDist'),   optsNg.MaxSnapDist   = []; end      % 允许的硬贴近距离（样点）；空则不限制

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
    [smax,ix] = max(score);
    xi_star = xi_list(ix);
    Ng_hat  = xi_star - N_hat;

    % 若启用“偏好 2 的幂次”，在近似等优集合中优先选择幂次 Ng
    if optsNg.PreferPow2
        keep = (score >= optsNg.SoftKeepRatio * smax);
        xi_keep = xi_list(keep);
        Ng_keep = xi_keep - N_hat;
        % 在 keep 集合里找 Pow2List 交集
        pow2_cands = intersect(Ng_keep, optsNg.Pow2List(:));
        if ~isempty(pow2_cands)
            % 从 keep∩pow2 中选分数最高者
            [~,best_idx] = max(score(keep & ismember(xi_list, N_hat + pow2_cands)));
            % best_idx 是 keep 掩码内的相对索引，转换到全局
            global_idx = find(keep);
            ix2 = global_idx(best_idx);
            xi_star = xi_list(ix2);
            Ng_hat  = xi_star - N_hat;
        else
            % 若 keep 集合内没有幂次，则在全体中选 Pow2List 中分数最高者
            xi_pow2 = N_hat + intersect((xi_list - N_hat), optsNg.Pow2List(:));
            [found, loc] = ismember(xi_pow2, xi_list);
            if any(found)
                sc_pow2 = score(loc(found));
                [~, j] = max(sc_pow2);
                xi_star = xi_pow2(find(found,1,'first') + j - 1);
                Ng_hat  = xi_star - N_hat;
            elseif ~isempty(optsNg.MaxSnapDist)
                % 或者，在全体中选与某个幂次最近且距离≤MaxSnapDist 的候选
                best_ix = ix; best_dist = inf;
                for v = optsNg.Pow2List(:).'
                    xi_t = N_hat + v;
                    [~,j] = min(abs(xi_list - xi_t));
                    dist = abs((xi_list(j) - N_hat) - v);
                    if dist <= optsNg.MaxSnapDist && score(j) > score(best_ix)
                        best_ix = j; best_dist = dist;
                    end
                end
                ix = best_ix; xi_star = xi_list(ix); Ng_hat = xi_star - N_hat;
            end
        end
    end

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
