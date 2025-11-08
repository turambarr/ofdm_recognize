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
    % Ng 偏好设置（可选：优先 2 的幂次，如 32）
    if ~isfield(opts,'NgPreferPow2'), opts.NgPreferPow2 = true; end
    if ~isfield(opts,'NgPow2List'),   opts.NgPow2List   = [16 32 64 128 256 512]; end
    if ~isfield(opts,'NgSoftKeepRatio'), opts.NgSoftKeepRatio = 0.98; end
    if ~isfield(opts,'NgMaxSnapDist'),   opts.NgMaxSnapDist   = []; end  % 为空不启用硬贴近

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
    % 传递 Ng 偏好选项
    ng_opts = struct('PreferPow2', opts.NgPreferPow2, ...
             'Pow2List',   opts.NgPow2List, ...
             'SoftKeepRatio', opts.NgSoftKeepRatio, ...
             'MaxSnapDist',   opts.NgMaxSnapDist);
    [Ng_hat, diagG] = paper_estimate_Ng(y_rs, N_hat, Fs_hat, opts.Np, opts.Td_ns, opts.Verbose, opts.NgTopK, ng_opts);
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
    R.all_Fs_candidates = all_results;  % 保存所有 Fs_guess 候选的结果
    R.best_candidate_idx = best_idx;    % 最佳候选索引
end

% —— 小工具：本文件局部使用（不生成新 .m 文件） ——
function R = fail_out(prefix, msg)
    fprintf(2,'[ERR] %s：%s\n', prefix, msg);
    R = struct('status','fail','errorMsg',[prefix '：' msg]);
end
