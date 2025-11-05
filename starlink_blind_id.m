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
    if ~isfield(opts,'Fs_guess'),opts.Fs_guess= 240e6; end         % 仅用于构造 Sr
    if ~isfield(opts,'p'),       opts.p       = 0.12; end          % Sr 相对宽度(±p)
    if ~isfield(opts,'nu_dB'),   opts.nu_dB   = 10; end            % 验证门限
    if ~isfield(opts,'M_N'),     opts.M_N     = []; end            % R0 使用样本数

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

        fprintf('[IO] 读到 %d 样点；format=%s\n', numel(y), meta.detected_format);
    catch ME
        R = fail_out('读取失败', ME.message); return;
    end

    % ---------- 2) 估 N、Nr、Fs ----------
    try
        [N_hat, Nr_hat, Fs_hat, diagN] = paper_estimate_N_and_Fs( ...
            y, Fr, opts.S_N, opts.Fs_guess, opts.p, opts.nu_dB, opts.M_N, opts.Verbose);
        fprintf('[OUT] N=%d, Nr=%d, Fs_hat=%.6f MHz\n', N_hat, Nr_hat, Fs_hat/1e6);
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
end

% —— 小工具：本文件局部使用（不生成新 .m 文件） ——
function R = fail_out(prefix, msg)
    fprintf(2,'[ERR] %s：%s\n', prefix, msg);
    R = struct('status','fail','errorMsg',[prefix '：' msg]);
end
