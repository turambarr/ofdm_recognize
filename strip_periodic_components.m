function [y_clean, info] = strip_periodic_components(y, opts)
% STRIP_PERIODIC_COMPONENTS  剥除信号中的未知导频/同步序列等“稳定重复”成分
%
% 该函数独立于 N/Fs 估计，纯数据驱动：
% - 先用“固定窗口+归一化”的线性自相关在一段 tau 范围内寻找显著周期；
% - 对每个显著 tau，按周期折叠/平均得到模板，再最小二乘缩放后从原信号中减去；
% - 可迭代移除多个显著周期，抑制干扰 N 估计的重复训练/同步结构；
%
% 用法:
%   [y2, info] = strip_periodic_components(y, opts)
%
% 输入:
%   y    : 复基带列向量
%   opts : 结构体（均为可选项）
%       .tau_min        最小滞后（样点），默认 16
%       .tau_max        最大滞后（样点），默认 min(16384, floor(numel(y)/4))
%       .maxPeaks       最多剥除的周期个数，默认 3
%       .thr_rho        归一化相关阈值（峰值），默认 0.2
%       .Kseg           多段一致性分段数（>=1，越大越稳），默认 3
%       .Mfrac          相关时使用的窗口长度 m0 = Mfrac * (有效长度)，默认 0.5
%       .recomputeEach  每剥除一个周期后是否重新计算 R(τ) 再找下一个，默认 true
%       .verbose        打印过程，默认 true
%
% 输出:
%   y_clean : 剥除重复成分后的信号
%   info    : 诊断结构体
%       .taus_detected      被剥除的周期列表
%       .rho_curve          归一化相关曲线（初始）
%       .taus_grid          计算 rho 的滞后网格
%       .removed_energy_pct 每次剥除的能量占比
%       .templates          每个 tau 的模板（cell）
%       .scales             每次拟合的复缩放因子 alpha
%       .y_len              样本长度
%
% 备注:
% - 该方法倾向剥除“在全局上反复出现且相位稳定”的成分，
%   比如重复的同步序列/导频模式；对非稳定/非周期数据影响较小。
% - 不需要任何先验参数（N/Fs/Ng 不参与）。
% - 若信号存在明显色噪/窄带，建议在外部先做预白化（可参考你的主流程）。
%
% 作者: 自动生成（独立模块）

    if nargin<2, opts = struct(); end
    y = y(:);
    L = numel(y);

    % 默认参数
    if ~isfield(opts,'tau_min'),       opts.tau_min = 16; end
    if ~isfield(opts,'tau_max'),       opts.tau_max = min(16384, floor(L/4)); end
    if ~isfield(opts,'maxPeaks'),      opts.maxPeaks = 3; end
    if ~isfield(opts,'thr_rho'),       opts.thr_rho = 0.20; end
    if ~isfield(opts,'Kseg'),          opts.Kseg = 3; end
    if ~isfield(opts,'Mfrac'),         opts.Mfrac = 0.5; end
    if ~isfield(opts,'recomputeEach'), opts.recomputeEach = true; end
    if ~isfield(opts,'verbose'),       opts.verbose = true; end

    y_work = y;
    removed_energy_pct = [];
    taus_removed = [];
    templates = {};
    scales = [];

    % 预先计算一次 rho 曲线（首次）
    [rho, taus] = rho_curve_multiseg(y_work, opts.tau_min, opts.tau_max, opts.Mfrac, opts.Kseg);

    if opts.verbose
        % 初始峰表（便于外部对比）
        [pks0, locs0] = findpeaks(rho, 'SortStr','descend');
        showK = min(8, numel(pks0));
        if showK>0
            fprintf('[strip] 初始前%d个周期峰: ', showK);
            for k=1:showK
                fprintf('#%d tau=%d rho=%.3f  ', k, taus(locs0(k)), pks0(k));
            end
            fprintf('\n');
        end
    end

    for it = 1:max(1, opts.maxPeaks)
        % 寻找峰值
        [pks, locs] = findpeaks(rho, 'MinPeakHeight', opts.thr_rho);
        if isempty(pks)
            if opts.verbose, fprintf('[strip] 没有高于阈值的周期峰，停止。
'); end
            break;
        end
        % 选择当前最强峰（也可考虑排除过近的重复峰）
        [~,ix] = max(pks);
        tau_sel = taus(locs(ix));
        if opts.verbose
            fprintf('[strip] 迭代 #%d: 选择 tau=%d (rho=%.3f) 进行剥除\n', it, tau_sel, pks(ix));
        end

        % 构造周期模板并估计缩放系数
        tmpl = estimate_period_template(y_work, tau_sel);
        if all(tmpl==0)
            if opts.verbose, fprintf('[strip] 模板全零，跳过。\n'); end
            break;
        end
        Yp = tile_template(tmpl, L);
        % 最小二乘尺度（复数）：argmin || y - alpha*Yp ||
        alpha = (Yp' * y_work) / max((Yp' * Yp), eps);
        if opts.verbose
            fprintf('[strip] tau=%d: 拟合缩放 alpha=%.4g%+.4gj\n', tau_sel, real(alpha), imag(alpha));
        end
        comp = alpha * Yp;

        % 能量统计
        E_before = sum(abs(y_work).^2) + eps;
        y_work2 = y_work - comp;
        E_after  = sum(abs(y_work2).^2) + eps;
        pct = max(0, min(1, (E_before - E_after)/E_before));
        if opts.verbose
            fprintf('[strip] tau=%d: 去除能量占比 ≈ %.2f%%\n', tau_sel, pct*100);
        end

        % 记录
        taus_removed(end+1,1) = tau_sel; %#ok<AGROW>
        removed_energy_pct(end+1,1) = pct; %#ok<AGROW>
        templates{end+1,1} = tmpl; %#ok<AGROW>
        scales(end+1,1) = alpha; %#ok<AGROW>
        y_work = y_work2;

        % 需要的话，重新计算 rho，继续下一轮
        if opts.recomputeEach
            [rho, taus] = rho_curve_multiseg(y_work, opts.tau_min, opts.tau_max, opts.Mfrac, opts.Kseg);
            if opts.verbose
                [pksR, locsR] = findpeaks(rho, 'SortStr','descend');
                showK2 = min(5, numel(pksR));
                if showK2>0
                    fprintf('[strip] 迭代后前%d峰: ', showK2);
                    for k2=1:showK2
                        fprintf('#%d tau=%d rho=%.3f  ', k2, taus(locsR(k2)), pksR(k2));
                    end
                    fprintf('\n');
                end
            end
        else
            % 否则从已有曲线中删除附近的峰，再看是否还有可剥除的
            mask_near = abs(taus - tau_sel) <= round(0.05 * tau_sel);
            rho(mask_near) = 0;
        end

        % 可选的停止条件：剥除能量很小则提前停
        if pct < 0.5/100  % <0.5%
            if opts.verbose, fprintf('[strip] 单次贡献 <0.5%%，提前停止。\n'); end
            break;
        end
    end

    y_clean = y_work;

    % 填充信息
    info = struct();
    info.taus_detected = taus_removed;
    info.removed_energy_pct = removed_energy_pct;
    info.templates = templates;
    info.scales = scales;
    info.y_len = L;
    info.rho_curve = rho;
    info.taus_grid = taus;
end

%% ================== 辅助函数 ==================
function [rho, taus] = rho_curve_multiseg(y, tmin, tmax, Mfrac, K)
% 归一化线性自相关（固定窗口）+ 多段一致性平均
    y = y(:); L = numel(y);
    tmin = max(1, round(tmin));
    tmax = min(L-2, round(tmax));
    taus = tmin:tmax;
    if isempty(taus), rho = zeros(0,1); return; end

    K = max(1, round(K));
    seg_len = floor(L / K);
    acc = 0; cnt = 0;
    for k=1:K
        s1 = (k-1)*seg_len + 1;
        if k<K, s2 = k*seg_len; else, s2 = L; end
        if s2 - s1 + 1 <= max(taus)+2, continue; end
        rho_k = rho_curve_single(y(s1:s2), taus, Mfrac);
        acc = acc + rho_k; cnt = cnt + 1;
    end
    if cnt>0
        rho = acc / cnt;
    else
        rho = rho_curve_single(y, taus, Mfrac);
    end
end

function rho = rho_curve_single(y, taus, Mfrac)
% 单段：固定 m0 + 归一化相关
    L = numel(y); taus = taus(:).';
    m0 = max(8, min([floor(Mfrac*(L-1)), L-1 - max(taus), L-1]));
    if m0 <= 0
        rho = zeros(numel(taus),1);
        return;
    end
    seg1 = y(1:m0);
    E1 = sum(abs(seg1).^2);
    rho = zeros(numel(taus),1);
    for i=1:numel(taus)
        t = taus(i); if t>=L-1, rho(i)=0; continue; end
        seg2 = y(1+t : t+m0);
        num = abs(sum(seg2 .* conj(seg1)));
        E2  = sum(abs(seg2).^2);
        den = sqrt(max(E1*E2, eps));
        rho(i) = num / den;
    end
end

function tmpl = estimate_period_template(y, tau)
% 按周期折叠平均得到模板（长度 = tau）
    y = y(:); L = numel(y);
    tau = round(tau);
    if tau<2 || tau>=L, tmpl = zeros(tau,1); return; end
    P = floor(L / tau);
    if P < 2, tmpl = zeros(tau,1); return; end
    % 构造矩阵: 每列为一个周期
    Y = reshape(y(1:tau*P), tau, P);
    tmpl = mean(Y, 2);
end

function yrep = tile_template(tmpl, L)
% 将模板重复铺满长度 L
    tau = numel(tmpl);
    reps = ceil(L / tau);
    yrep = repmat(tmpl, reps, 1);
    yrep = yrep(1:L);
end
