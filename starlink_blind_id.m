function R = starlink_blind_id(inputFile, Fr, opts)
% STARLINK_BLIND_ID  —— 极简主流程（只用 3 个函数文件）
% 1) read_iq_autodetect 读取交错IQ（自动识别 int16/float32）
% 2) paper_estimate_N_and_Fs 估计 N、Nr、Fs
% 3) （内联 interp1 重采样到 Fs_hat）
% 4) paper_estimate_Ng 估计 Ng
%
% 仅识别 N、Fs、Ng（及派生的 Tsym、F 等），不做其它步骤。
%
% 用法:
%   Fr   = 409.6e6;   % 采样率（Hz），用你采集链路的真实 Fr
%   opts = struct();  % 可留空使用默认
%   R = starlink_blind_id('IQ_xxx.dat', Fr, opts);

    if nargin<2
        fprintf('用法: R = starlink_blind_id(inputFile, Fr, opts)\n');
        R = struct('status','fail','errorMsg','参数不足'); return;
    end
    if nargin<3 || isempty(opts), opts = struct(); end

    % ---------- 预设参数（可被 opts 覆盖） ----------
    if ~isfield(opts,'format'),    opts.format = 'auto'; end   % 'auto'|'int16_iq'|'float32_iq'
    if ~isfield(opts,'scale'),     opts.scale  = 1/32768; end
    if ~isfield(opts,'endian'),    opts.endian = 'ieee-le'; end
    if ~isfield(opts,'order'),     opts.order  = 'IQ'; end
    if ~isfield(opts,'MaxSamples'),opts.MaxSamples = []; end   % 为空则全用
    if ~isfield(opts,'Detrend'),   opts.Detrend = true; end    % 去均值
    if ~isfield(opts,'Verbose'),   opts.Verbose = false; end   % 显示候选过程
    if ~isfield(opts,'NgTopK'),    opts.NgTopK = 10; end       % 显示前 K 个 Ng 候选
    % 周期性成分剥除（独立预处理）
    if ~isfield(opts,'EnableStripPeriodic'), opts.EnableStripPeriodic = false; end
    if ~isfield(opts,'StripPeriodicOpts'),   opts.StripPeriodicOpts   = struct(); end
    % 调试输出
    if ~isfield(opts,'DebugLevel'),  opts.DebugLevel = 0; end   % 0:无 1:关键 2:详细
    if ~isfield(opts,'DebugTopK'),   opts.DebugTopK  = 8; end

    % 论文参数（N/Fs）                
    if ~isfield(opts,'S_N'),     opts.S_N     = 2.^(8:13); end     % N 候选
    if ~isfield(opts,'Fs_guess'),opts.Fs_guess= []; end         % 仅用于构造 Sr（可为数组）
    if ~isfield(opts,'p'),       opts.p       = 0.05; end          % Sr 相对宽度(±p)
    if ~isfield(opts,'nu_dB'),   opts.nu_dB   = 10; end            % 验证门限
    if ~isfield(opts,'M_N'),     opts.M_N     = []; end            % R0 使用样本数
    
    % CP 成对峰验证参数（新增）
    if ~isfield(opts,'EnableCPCheck'),  opts.EnableCPCheck = true; end   % 启用 CP 检查
    if ~isfield(opts,'CPRatioMin'),     opts.CPRatioMin = 0.15; end      % CP 副峰最小比例
    if ~isfield(opts,'CPTolerancePct'), opts.CPTolerancePct = 0.15; end  % CP 位置容差
    if ~isfield(opts,'CPPenaltyDB'),    opts.CPPenaltyDB = 15; end       % 无 CP 峰惩罚 dB
    
    % 确保 Fs_guess 是行向量
    opts.Fs_guess = opts.Fs_guess(:).';
    num_Fs_candidates = numel(opts.Fs_guess);

    % 论文参数（Ng）
    if ~isfield(opts,'Np'),      opts.Np      = 4; end
    if ~isfield(opts,'Td_ns'),   opts.Td_ns   = 108; end

    % ---------- 1) 读取 IQ ----------
    try
        [y, meta] = read_iq_autodetect(inputFile, opts);
        if ~isempty(opts.MaxSamples)
            y = y(1:min(end, floor(opts.MaxSamples)));
        end

        if opts.Detrend
            y = y - mean(y);
        end
        % RMS 归一化
        r = sqrt(mean(abs(y).^2)); if r>0, y = y / r; end

        % ===== 频谱预白化 =====
        Yf = fft(y);
        Px = abs(Yf).^2;
        % 101点移动平均平滑功率谱
        Px_smooth = movmean(Px, 101);
        W = 1 ./ sqrt(Px_smooth + 1e-12);
        % 限定白化滤波器动态范围，避免极端放大
        W = min(max(W, 0.3), 3);
        Yw = Yf .* W;
        y_white = ifft(Yw);
        y = y_white;

        % —— Debug：剥除前的周期峰列表（快速RHO）——
        if opts.DebugLevel>=1
            try
                Kseg_dbg = 3; Mfrac_dbg = 0.5; tmin_dbg = 16; tmax_dbg = min(20000, floor(numel(y)/4));
                [rho0, taus0] = local_rho_curve_multiseg(y, tmin_dbg, tmax_dbg, Mfrac_dbg, Kseg_dbg);
                [pks, locs] = findpeaks(rho0);
                [pks_sorted, ord] = sort(pks, 'descend');
                Ltop = min(opts.DebugTopK, numel(ord));
                fprintf('[DBG] RHO-PEAKS(before strip): ');
                for k=1:Ltop
                    fprintf('#%d tau=%d rho=%.3f  ', k, taus0(locs(ord(k))), pks_sorted(k));
                end
                fprintf('\n');
            catch
            end
        end

        % —— 可选：剥除未知导频/同步等稳定重复成分 ——
        strip_info = struct();
        if isfield(opts,'EnableStripPeriodic') && opts.EnableStripPeriodic
            sp = opts.StripPeriodicOpts;
            if ~isfield(sp,'thr_rho'),   sp.thr_rho = 0.22; end     % 显著周期阈值（保守→0.25；激进→0.18）
            if ~isfield(sp,'Kseg'),      sp.Kseg    = 5;    end     % 多段一致性，增强稳健性
            if ~isfield(sp,'maxPeaks'),  sp.maxPeaks= 2;    end     % 剥除最多2类强重复
            if ~isfield(sp,'Mfrac'),     sp.Mfrac   = 0.5;  end     % 相关窗口相对长度
            if ~isfield(sp,'recomputeEach'), sp.recomputeEach = true; end
            if ~isfield(sp,'verbose'),   sp.verbose = logical(opts.Verbose); end
            [y_stripped, strip_info] = strip_periodic_components(y, sp);
            if ~isempty(strip_info)
                fprintf('[STRIP] 剥除周期数=%d，总能量去除≈%.2f%%，taus=%s\n', ...
                    numel(strip_info.taus_detected), 100*sum(strip_info.removed_energy_pct), mat2str(strip_info.taus_detected.'));
            end
            y = y_stripped;
            % —— Debug：剥除后的周期峰列表 ——
            if opts.DebugLevel>=1
                try
                    Kseg_dbg = 3; Mfrac_dbg = 0.5; tmin_dbg = 16; tmax_dbg = min(20000, floor(numel(y)/4));
                    [rho1, taus1] = local_rho_curve_multiseg(y, tmin_dbg, tmax_dbg, Mfrac_dbg, Kseg_dbg);
                    [pks, locs] = findpeaks(rho1);
                    [pks_sorted, ord] = sort(pks, 'descend');
                    Ltop = min(opts.DebugTopK, numel(ord));
                    fprintf('[DBG] RHO-PEAKS(after strip):  ');
                    for k=1:Ltop
                        fprintf('#%d tau=%d rho=%.3f  ', k, taus1(locs(ord(k))), pks_sorted(k));
                    end
                    fprintf('\n');
                catch
                end
            end
        end

        fprintf('[IO] 读到 %d 样点；format=%s\n', numel(y), meta.detected_format);
    catch ME
        R = fail_out('读取失败', ME.message); return;
    end

    % ---------- 2) 估 N、Nr、Fs（支持多个 Fs_guess 候选） ----------
    try
        fprintf('\n========== 开始估计 N 和 Fs（共 %d 个 Fs_guess 候选）==========\n', num_Fs_candidates);
        
        % 为每个 Fs_guess 候选值进行估计
        all_results = cell(num_Fs_candidates, 1);
        for fs_idx = 1:num_Fs_candidates
            Fs_guess_current = opts.Fs_guess(fs_idx);
            fprintf('\n---------- Fs_guess 候选 #%d: %.3f MHz ----------\n', fs_idx, Fs_guess_current/1e6);
            
            % 构建 CP 检查参数
            cp_opts = struct('EnableCPCheck', opts.EnableCPCheck, ...
                             'CPRatioMin', opts.CPRatioMin, ...
                             'CPTolerancePct', opts.CPTolerancePct, ...
                             'CPPenaltyDB', opts.CPPenaltyDB);
            
            [N_hat, Nr_hat, Fs_hat, diagN] = paper_estimate_N_and_Fs( ...
                y, Fr, opts.S_N, Fs_guess_current, opts.p, opts.nu_dB, opts.M_N, opts.Verbose, cp_opts);
            
            fprintf('[OUT #%d] N=%d, Nr=%d, Fs_hat=%.6f MHz (基于 Fs_guess=%.3f MHz)\n', ...
                fs_idx, N_hat, Nr_hat, Fs_hat/1e6, Fs_guess_current/1e6);
            
            % 保存当前候选的结果
            all_results{fs_idx} = struct('Fs_guess', Fs_guess_current, ...
                                          'N_hat', N_hat, ...
                                          'Nr_hat', Nr_hat, ...
                                          'Fs_hat', Fs_hat, ...
                                          'diagN', diagN);
        end
        
        fprintf('\n========== 所有 Fs_guess 候选结果汇总 ==========\n');
        for fs_idx = 1:num_Fs_candidates
            res = all_results{fs_idx};
            cp_info = '';
            if isfield(res.diagN, 'cp_ratio_best')
                cp_info = sprintf(', CP=%.1f%% (contrast=%.1f%%)', res.diagN.cp_ratio_best*100, 100*getfield(res.diagN, 'cp_contrast_best'));
            end
            fprintf('[候选#%d] Fs_guess=%.3f MHz → N=%d, Nr=%d, Fs_hat=%.6f MHz, dyn=%.1f dB%s\n', ...
                fs_idx, res.Fs_guess/1e6, res.N_hat, res.Nr_hat, res.Fs_hat/1e6, ...
                res.diagN.dyn_dB_best, cp_info);
        end
        
        % 选择最佳结果（如果启用 CP，基于 cp_score；否则基于 dyn_dB）
        score_values = zeros(num_Fs_candidates, 1);
        for fs_idx = 1:num_Fs_candidates
            diagN_fs = all_results{fs_idx}.diagN;
            if opts.EnableCPCheck && isfield(diagN_fs, 'cp_score')
                % 使用 CP 调整后的得分（选最高的那个 N 候选的 cp_score）
                [~, best_N_idx] = max(diagN_fs.cp_score);
                score_values(fs_idx) = diagN_fs.cp_score(best_N_idx);
            else
                score_values(fs_idx) = diagN_fs.dyn_dB_best;
            end
        end
        [best_score, best_idx] = max(score_values);
        
        if opts.EnableCPCheck
            fprintf('\n[最佳] 选择候选 #%d (Fs_guess=%.3f MHz, CP验证得分=%.1f dB 最高)\n', ...
                best_idx, opts.Fs_guess(best_idx)/1e6, best_score);
        else
            fprintf('\n[最佳] 选择候选 #%d (Fs_guess=%.3f MHz, dyn=%.1f dB 最高)\n', ...
                best_idx, opts.Fs_guess(best_idx)/1e6, best_score);
        end
        
        % 使用最佳结果继续后续处理
        N_hat  = all_results{best_idx}.N_hat;
        Nr_hat = all_results{best_idx}.Nr_hat;
        Fs_hat = all_results{best_idx}.Fs_hat;
        diagN  = all_results{best_idx}.diagN;
        
    catch ME
        R = fail_out('[N/Fs] 估计失败', ME.message); return;
    end

    % ---------- 3) 重采样到 Fs_hat（仅用内置 interp1，不引入新文件） ----------
    try
        if abs(Fs_hat - Fr) < 1e-9*Fr
            y_rs = y;
        else
            L  = numel(y);
            t1 = (0:L-1).' / Fr;
            L2 = max(1, round(L * Fs_hat / Fr));
            t2 = (0:L2-1).' / Fs_hat;
            y_rs = interp1(t1, y, t2, 'linear', 'extrap');  % 线性插值够用
        end
    catch ME
        R = fail_out('重采样失败', ME.message); return;
    end

    % ---------- 4) 估 Ng ----------
    try
        [Ng_hat, diagG] = paper_estimate_Ng(y_rs, N_hat, Fs_hat, opts.Np, opts.Td_ns, opts.Verbose, opts.NgTopK);
        Tsamp_sym = N_hat + Ng_hat;
        Tsym_hat  = Tsamp_sym / Fs_hat;
        F_sub     = Fs_hat / N_hat;
        fprintf('[OUT] Ng=%d → Tsym≈%.3f us, F≈%.3f kHz\n', Ng_hat, Tsym_hat*1e6, F_sub/1e3);
    catch ME
        R = fail_out('[Ng] 估计失败', ME.message); return;
    end

    % ---------- 汇总输出 ----------
    R = struct();
    R.status = 'ok'; R.errorMsg = '';
    R.meta = meta;
    R.Fr = Fr; R.Fs_hat = Fs_hat;
    R.N_hat = N_hat; R.Nr_hat = Nr_hat;
    R.Ng_hat = Ng_hat;
    R.Tsym_samples = N_hat + Ng_hat;
    R.Tsym = (N_hat + Ng_hat)/Fs_hat;
    R.F = Fs_hat / N_hat;      % 子载波间隔
    R.diagN = diagN; R.diagG = diagG;
    if exist('strip_info','var') && ~isempty(strip_info)
        R.strip_info = strip_info;
    else
        R.strip_info = [];
    end
    R.all_Fs_candidates = all_results;  % 保存所有 Fs_guess 候选的结果
    R.best_candidate_idx = best_idx;    % 最佳候选索引
end

% —— 小工具：本文件局部使用（不生成新 .m 文件） ——
function R = fail_out(prefix, msg)
    fprintf(2,'[ERR] %s：%s\n', prefix, msg);
    R = struct('status','fail','errorMsg',[prefix '：' msg]);
end

function [rho, taus] = local_rho_curve_multiseg(y, tmin, tmax, Mfrac, K)
% 供调试打印使用：固定窗口+归一化相关 + 多段平均
    y = y(:); L = numel(y);
    tmin = max(1, round(tmin)); tmax = min(L-2, round(tmax));
    taus = tmin:tmax;
    if isempty(taus), rho = zeros(0,1); return; end
    K = max(1, round(K));
    seg_len = floor(L / K);
    acc = 0; cnt=0;
    for k=1:K
        s1 = (k-1)*seg_len + 1; if k<K, s2 = k*seg_len; else, s2 = L; end
        if s2 - s1 + 1 <= max(taus)+2, continue; end
        acc = acc + local_rho_curve_single(y(s1:s2), taus, Mfrac);
        cnt = cnt + 1;
    end
    if cnt>0, rho = acc/cnt; else, rho = local_rho_curve_single(y, taus, Mfrac); end
end

function rho = local_rho_curve_single(y, taus, Mfrac)
    L = numel(y); taus = taus(:).';
    m0 = max(8, min([floor(Mfrac*(L-1)), L-1 - max(taus), L-1]));
    if m0<=0, rho = zeros(numel(taus),1); return; end
    seg1 = y(1:m0); E1 = sum(abs(seg1).^2);
    rho = zeros(numel(taus),1);
    for i=1:numel(taus)
        t = taus(i); if t>=L-1, rho(i)=0; continue; end
        seg2 = y(1+t:t+m0);
        num = abs(sum(seg2 .* conj(seg1)));
        E2  = sum(abs(seg2).^2);
        rho(i) = num / max(sqrt(E1*E2), eps);
    end
end
