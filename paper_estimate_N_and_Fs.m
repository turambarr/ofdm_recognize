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
R0mag = r0_curve_abs(y, taus, M);  % |R0(τ)|（式(11)）

% 如果启用 CP 检查，预先计算全局 R0 曲线（扩展范围以包含 CP 副峰位置）
if opts.EnableCPCheck
    % 扩展 taus 范围，确保能覆盖可能的 Ng 位置
    tau_min = 1;
    tau_max = max(taus);
    taus_extended = unique([tau_min:tau_max, taus(:)']);
    R0mag_extended = r0_curve_abs(y, taus_extended, M);
else
    taus_extended = taus;
    R0mag_extended = R0mag;
end

% 分桶找峰 + 每候选的 dyn(dB) + 改进 CP 验证
peak_vals = -inf(numel(S_N),1);
tau_arg   = zeros(numel(S_N),1);
dyn_vec   = -inf(numel(S_N),1);
cp_ratio  = zeros(numel(S_N),1);        % |R(Ng)| / |R(N+Ng)| （截断≤1）
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
    dyn_vec(k) = 20*log10(max(seg) / max(min(seg), eps));
    
    % ===== 改进 CP 成对峰验证：只检查 Ng = tau_arg - N_cand =====
    if opts.EnableCPCheck
        Ng_cand = tau_arg(k) - N_cand;    % 由峰值推导的 CP 长度候选
        Ng_cand_vec(k) = Ng_cand;
        if Ng_cand > 0
            % 合理 Ng 比例判断
            frac = Ng_cand / N_cand;
            ng_range_ok(k) = (frac >= opts.NgFracMin) && (frac <= opts.NgFracMax);
            % 在扩展自相关中直接取该 Ng 的相关值（若存在）
            if Ng_cand <= max(taus_extended)
                cp_peak = R0mag_extended( taus_extended == Ng_cand );
                if isempty(cp_peak), cp_peak = 0; end
            else
                cp_peak = 0;
            end
            cp_ratio(k) = min( cp_peak / max(mx, eps), 1.0 ); % 截断到 ≤1
            cp_pass_flag(k) = (cp_ratio(k) >= opts.CPRatioMin) && ng_range_ok(k);
            % 得分策略：通过则 dyn + 提升；不通过则 dyn - 惩罚
            if cp_pass_flag(k)
                cp_score(k) = dyn_vec(k) + opts.CPBoostWeight * cp_ratio(k); % 增强分
            else
                penalty = opts.CPPenaltyDB;
                if ~ng_range_ok(k), penalty = penalty + opts.NgFracPenaltyDB; end
                cp_score(k) = dyn_vec(k) - penalty * (1 - cp_ratio(k));
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

% 验证门限：使用 CP 验证得分选择最佳候选
% 优先选择 cp_score 最高且通过 nu_dB 门限的候选
cp_score_f = cp_score;
while true
    [~, i_best] = max(cp_score_f);
    if ~isfinite(cp_score_f(i_best))
        error('没有通过门限的候选；请增大 p 或 M，或放宽 CP 检查参数。');
    end
    if dyn_vec(i_best) >= nu_dB
        dyn_dB_best = dyn_vec(i_best);
        cp_ratio_best = cp_ratio(i_best);
        break;
    else
        cp_score_f(i_best) = -inf;
    end
end

N_hat  = S_N(i_best);
Nr_hat = tau_arg(i_best);
Ng_best = max(Nr_hat - N_hat, 0);  % 推导的 Ng（仅用于诊断；真正 Ng 在后续 Ng 估计步骤得）
Fs_hat = round( N_hat * (Fr / max(Nr_hat,1)) ); % 仍保留原公式（假设 Nr≈N+Ng）

if Verbose
    if opts.EnableCPCheck
        fprintf('[N/Fs] Chosen: N=%d, Nr=%d, (Ng≈%d, Ng/N=%.2f%%) Fs_hat=%.6f MHz dyn=%.1f dB cp_ratio=%.1f%% score=%.1f\n', ...
            N_hat, Nr_hat, Ng_best, 100*Ng_best/max(N_hat,1), Fs_hat/1e6, dyn_dB_best, cp_ratio_best*100, cp_score(i_best));
    else
        fprintf('[N/Fs] Chosen: N=%d, Nr=%d, Fs_hat=%.6f MHz, dyn=%.1f dB\n', ...
            N_hat, Nr_hat, Fs_hat/1e6, dyn_dB_best);
    end
end

diagN = struct('taus',taus,'R0mag',R0mag,'S',S_N,'peak_vals',peak_vals, ...
    'tau_arg',tau_arg,'dyn_dB',dyn_vec,'dyn_dB_best',dyn_dB_best, ...
    'cp_ratio',cp_ratio,'cp_score',cp_score,'cp_ratio_best',cp_ratio_best, ...
    'Ng_cand',Ng_cand_vec,'cp_pass',cp_pass_flag,'Ng_range_ok',ng_range_ok, ...
    'taus_extended',taus_extended,'R0mag_extended',R0mag_extended);
end

function s = ternary_str(c,a,b)
if c, s=a; else, s=b; end
end

function mag = r0_curve_abs(y, taus, M)
L = numel(y); mag = zeros(size(taus));
for i=1:numel(taus)
    t=taus(i); if t>=L-1, mag(i)=0; continue; end
    m=min(M, L-1-t);
    seg = y(1+t:m+t) .* conj(y(1:m));
    mag(i) = abs(sum(seg))/m;
end
end
