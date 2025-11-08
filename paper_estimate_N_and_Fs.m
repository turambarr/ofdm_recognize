function [N_hat, Nr_hat, Fs_hat, diagN] = paper_estimate_N_and_Fs(y, Fr, S_N, Fs_guess, p, nu_dB, M, Verbose, opts)
% 构造接收机域区间 Sr（式(17)-(19)），扫 |R0(τ)|（式(11)近似），做 ν(dB) 验证（式(20)），得 N、Nr、Fs（式(22)）
% 改进版：加入更严格的 CP 成对峰验证（基于 Ng = τ - N），排除虚假强周期。
y = y(:); L = numel(y);
if nargin<7 || isempty(M), M = min(L-1, max(round(0.015*Fr), round(0.5*L))); end
if nargin<8 || isempty(Verbose), Verbose = false; end
if nargin<9 || isempty(opts), opts = struct(); end

% CP 验证参数（新逻辑只检查 Ng = τ - N）
if ~isfield(opts,'EnableCPCheck'),    opts.EnableCPCheck    = true; end     % 是否启用 CP 检查
if ~isfield(opts,'CPRatioMin'),       opts.CPRatioMin       = 0.12; end     % 认为“存在 CP”所需的相关强度比阈值（≤1）
if ~isfield(opts,'CPPenaltyDB'),      opts.CPPenaltyDB      = 15;   end     % 缺失 CP 惩罚
if ~isfield(opts,'NgFracMin'),        opts.NgFracMin        = 1/64; end     % 合理 Ng/N 最小值
if ~isfield(opts,'NgFracMax'),        opts.NgFracMax        = 1/4;  end     % 合理 Ng/N 最大值
if ~isfield(opts,'NgFracPenaltyDB'),  opts.NgFracPenaltyDB  = 10;   end     % Ng 比例不合理额外惩罚
if ~isfield(opts,'CPBoostWeight'),    opts.CPBoostWeight    = 5;    end     % 通过 CP 时给 dyn 的加权提升（经验）
if ~isfield(opts,'Kseg'),             opts.Kseg             = 3;    end     % 多段一致性分段数（>=1）
if ~isfield(opts,'CPTolerancePct'),   opts.CPTolerancePct   = 0.08; end    % CP 窗口相对 N 的容差

eta = Fr / Fs_guess;

% 构造 Sr & 分桶
cells = cell(numel(S_N),1); Sr = [];
for k=1:numel(S_N)
    b=S_N(k);
    lo=max(1, floor(b*eta*(1-p)));
    hi=max(lo, ceil(b*eta*(1+p)));
    cells{k}=lo:hi; Sr=[Sr, lo:hi]; 
end
taus = unique(Sr);
R0mag = r0_curve_abs_multiseg(y, taus, M, max(1, round(opts.Kseg)));  % 归一化相关 + 多段平均

% 如果启用 CP 检查，预先计算全局 R0 曲线（扩展范围以包含 CP 副峰位置）
if opts.EnableCPCheck
    % 扩展 taus 范围，确保能覆盖可能的 Ng 位置
    tau_min = 1;
    tau_max = max(taus);
    taus_extended = unique([tau_min:tau_max, taus(:)']);
    R0mag_extended = r0_curve_abs_multiseg(y, taus_extended, M, max(1, round(opts.Kseg))); % 多段平均
else
    taus_extended = taus;
    R0mag_extended = R0mag;
end

function mag = r0_curve_abs_multiseg(y, taus, M, K)
% 将序列切分为 K 段，分别计算归一化相关并取平均，提升起止位置稳健性
L = numel(y);
K = max(1, round(K));
if K==1
    mag = r0_curve_abs(y, taus, M);
    return;
end
seg_len = floor(L / K);
acc = 0;
cnt = 0;
for k=1:K
    s1 = (k-1)*seg_len + 1;
    if k<K, s2 = k*seg_len; else, s2 = L; end
    if s2 - s1 + 1 <= max(taus)+2, continue; end
    mag_k = r0_curve_abs(y(s1:s2), taus, M);
    acc = acc + mag_k;
    cnt = cnt + 1;
end
if cnt>0
    mag = acc / cnt;
else
    mag = r0_curve_abs(y, taus, M);
end
end

% 分桶找峰 + 每候选的 dyn(dB) + 改进 CP 验证
peak_vals = -inf(numel(S_N),1);
tau_arg   = zeros(numel(S_N),1);
dyn_vec   = -inf(numel(S_N),1);
cp_ratio  = zeros(numel(S_N),1);        % 绝对比 |R(Ng)|/|R(N+Ng)|
cp_contrast_vec = zeros(numel(S_N),1);  % 对比度 (副峰-局部中位)/主峰
cp_score  = -inf(numel(S_N),1);         % 最终用于选择的得分
Ng_cand_vec = -ones(numel(S_N),1);      % 每个候选对应的 Ng (=tau_arg - N)
cp_pass_flag = false(numel(S_N),1);     % CP 通过标记
ng_range_ok  = false(numel(S_N),1);     % Ng 比例是否在合理范围

for k=1:numel(S_N)
    N_cand = S_N(k);
    tset = cells{k}; if isempty(tset), continue; end
    mask = ismember(taus, tset);
    seg  = R0mag(mask);
    [mx, ix] = max(seg);
    peak_vals(k) = mx; 
    tau_arg(k) = tset(ix);
    % 用中位数作底噪，稳健动态范围
    base = median(seg);
        % 稳健动态范围：主峰相对窗口中位数
        dyn_vec(k) = 20*log10(max(seg) / max(median(seg), eps));
    
    % ===== 改进 CP 成对峰验证：只检查 Ng = tau_arg - N_cand =====
    if opts.EnableCPCheck
        Ng_cand = tau_arg(k) - N_cand;    % 由峰值推导的 CP 长度候选
        Ng_cand_vec(k) = Ng_cand;
        if Ng_cand > 0
            % 合理 Ng 比例判断
            frac = Ng_cand / N_cand;
            ng_range_ok(k) = (frac >= opts.NgFracMin) && (frac <= opts.NgFracMax);
            % 在 Ng±tol 窗内取局部峰（tol 由 opts.CPTolerancePct 指定）
            tol = max(1, round(opts.CPTolerancePct * N_cand));
            lo_cp = max(1, Ng_cand - tol); hi_cp = min(max(taus_extended), Ng_cand + tol);
            mask_cp = (taus_extended>=lo_cp & taus_extended<=hi_cp);
            if any(mask_cp)
                win_vals = R0mag_extended(mask_cp);
                cp_peak = max(win_vals);
                cp_med  = median(win_vals);
            else
                cp_peak = 0; cp_med = 0;
            end
            cp_ratio(k) = min( cp_peak / max(mx, eps), 1.0 ); % 绝对比
            cp_contrast_vec(k) = max(cp_peak - cp_med, 0) / max(mx, eps); % 对比度
            cp_pass_flag(k) = (cp_contrast_vec(k) >= opts.CPRatioMin) && ng_range_ok(k);
            % 得分策略：通过则 dyn + 提升；不通过则 dyn - 惩罚
            if cp_pass_flag(k)
                cp_score(k) = dyn_vec(k) + opts.CPBoostWeight * cp_contrast_vec(k); % 用对比度增强分
            else
                penalty = opts.CPPenaltyDB;
                if ~ng_range_ok(k), penalty = penalty + opts.NgFracPenaltyDB; end
                cp_score(k) = dyn_vec(k) - penalty * (1 - min(cp_contrast_vec(k),1));
            end
        else
            % Ng <=0 不合理
            cp_ratio(k) = 0;
            cp_score(k) = dyn_vec(k) - (opts.CPPenaltyDB + opts.NgFracPenaltyDB);
        end
    else
        % 未启用 CP 检查：得分即 dyn
        cp_score(k) = dyn_vec(k);
    end
    
    if Verbose
        passflag = dyn_vec(k) >= nu_dB;
        cp_flag = '';
        if opts.EnableCPCheck
            if cp_pass_flag(k)
                cp_flag = sprintf(' [CP✓ ratio=%.1f%% Ng=%d (%.2f%%)]', cp_ratio(k)*100, Ng_cand_vec(k), 100*Ng_cand_vec(k)/max(N_cand,1));
            else
                cp_flag = sprintf(' [CP✗ ratio=%.1f%% Ng=%d (%.2f%%)%s]', cp_ratio(k)*100, Ng_cand_vec(k), 100*max(Ng_cand_vec(k),0)/max(N_cand,1), ternary_str(ng_range_ok(k),'',' RANGE!')); 
            end
        end
        fprintf('[N/Fs] N=%d: tau=[%d,%d], peak=%.4g@%d, dyn=%.1f dB score=%.1f%s%s\n', ...
            N_cand, tset(1), tset(end), mx, tau_arg(k), dyn_vec(k), cp_score(k), ternary_str(passflag,' PASS',''), cp_flag);
    end
end

% 验证门限与回退策略：
% 1) 优先在 dyn >= nu_dB 的集合内选 cp_score 最大者；
% 2) 若无一满足门限，则回退为“全体中 cp_score 最大者”（并提示采用回退路径）。
idx_pass_dyn = find(dyn_vec >= nu_dB);
if ~isempty(idx_pass_dyn)
    [~, rel] = max(cp_score(idx_pass_dyn));
    i_best = idx_pass_dyn(rel);
    used_fallback = false;
else
    [~, i_best] = max(cp_score);
    used_fallback = true;
end
dyn_dB_best = dyn_vec(i_best);
cp_ratio_best = cp_ratio(i_best);
cp_contrast_best = cp_contrast_vec(i_best);

N_hat  = S_N(i_best);
Nr_hat = tau_arg(i_best);
Ng_best = max(Nr_hat - N_hat, 0);  % 推导的 Ng（仅用于诊断；真正 Ng 在后续 Ng 估计步骤得）
Fs_hat = round( N_hat * (Fr / max(Nr_hat,1)) ); % 仍保留原公式（假设 Nr≈N+Ng）

if Verbose
    if opts.EnableCPCheck
        if exist('used_fallback','var') && used_fallback
            fbmsg = ' (fallback: no dyn>=nu)';
        else
            fbmsg = '';
        end
        fprintf('[N/Fs] Chosen: N=%d, Nr=%d, (Ng≈%d, Ng/N=%.2f%%) Fs_hat=%.6f MHz dyn=%.1f dB cp_ratio=%.1f%% cp_contrast=%.3f score=%.1f%s\n', ...
            N_hat, Nr_hat, Ng_best, 100*Ng_best/max(N_hat,1), Fs_hat/1e6, dyn_dB_best, cp_ratio_best*100, cp_contrast_best, cp_score(i_best), fbmsg);
    else
        fprintf('[N/Fs] Chosen: N=%d, Nr=%d, Fs_hat=%.6f MHz, dyn=%.1f dB\n', ...
            N_hat, Nr_hat, Fs_hat/1e6, dyn_dB_best);
    end
end

diagN = struct('taus',taus,'R0mag',R0mag,'S',S_N,'peak_vals',peak_vals, ...
    'tau_arg',tau_arg,'dyn_dB',dyn_vec,'dyn_dB_best',dyn_dB_best, ...
    'cp_ratio',cp_ratio,'cp_contrast',cp_contrast_vec,'cp_score',cp_score,'cp_ratio_best',cp_ratio_best,'cp_contrast_best',cp_contrast_best, ...
    'Ng_cand',Ng_cand_vec,'cp_pass',cp_pass_flag,'Ng_range_ok',ng_range_ok, ...
    'taus_extended',taus_extended,'R0mag_extended',R0mag_extended);
end

function s = ternary_str(c,a,b)
if c, s=a; else, s=b; end
end

function mag = r0_curve_abs(y, taus, M)
% 固定窗口长度+归一化相关
L = numel(y); mag = zeros(size(taus));
taus = taus(:).';
if isempty(taus), return; end
m0 = min([M, L-1-max(taus), L-1]);
if m0 <= 0, return; end
seg1 = y(1:m0);
E1 = sum(abs(seg1).^2);
for i=1:numel(taus)
    t = taus(i); if t>=L-1, mag(i)=0; continue; end
    seg2 = y(1+t : t+m0);
    num = abs(sum(seg2 .* conj(seg1)));
    E2  = sum(abs(seg2).^2);
    den = sqrt(max(E1*E2, eps));
    mag(i) = num / den;
end
end
