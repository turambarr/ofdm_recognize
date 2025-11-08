function estimate_Fs_guess()
clc; close all;

%% ===== 1) 配置 =====
Fr = 409.6e6;              % ← 必填：接收机采样率(Hz)，如 20e6。留空会弹窗让你输入
infile = 'test3.dat'; % 默认数据文件（int16 I/Q 交织）；不存在会弹窗
iq_order = 'iq';      % 'iq' 或 'qi'（I、Q 顺序）
eta_grid = [0.80 0.85 0.90 0.95];   % 常见占用率（可改）
S_N = [256 512 1024 1536 2048 3072 4096 6144 8192]; % N 候选（可加/减）
p = 0.20;             % τ 窗口相对宽度（±p），粗扫用 0.2 稳一点
nu_soft = 6;          % 动态范围“软阈”用于排序（dB，不硬裁)

%% ===== 2) 读文件 & 采样率 =====
if ~exist(infile,'file')
    [fn, fp] = uigetfile({'*.dat;*.bin;*.*','All data files (*.dat, *.bin, *.*)'}, ...
                          '请选择 IQ 数据文件（int16 I/Q 交织）');
    if isequal(fn,0), error('未选择文件'); end
    infile = fullfile(fp,fn);
end
raw = read_int16(infile);
if mod(numel(raw),2)~=0
    error('文件样本数不是偶数，无法按 I/Q 交织读取。');
end
switch lower(iq_order)
    case 'iq', I = double(raw(1:2:end)); Q = double(raw(2:2:end));
    case 'qi', Q = double(raw(1:2:end)); I = double(raw(2:2:end));
    otherwise, error('iq_order 仅支持 ''iq'' 或 ''qi''');
end
y = complex(I,Q); clear I Q raw;
y = y(:) - mean(y);
L = numel(y);
fprintf('读取到 %.3f M 复样点。\n', L/1e6);

if isempty(Fr)
    prompt = {'请输入接收机采样率 Fr（Hz）:'};
    answ = inputdlg(prompt,'Fr 必填',1,{'20e6'});
    if isempty(answ), error('未填写 Fr。'); end
    Fr = str2double(answ{1});
    if ~isfinite(Fr) || Fr<=0, error('Fr 非法。'); end
end

%% ===== 3) 估占用带宽 B_occ（两法并行，取稳健融合）=====
% periodogram（中心化）
nfft = max(2^12, 2^nextpow2(min(L, 2^20)));
[Px, f] = periodogram(y, [], nfft, Fr, 'centered');   % 用默认窗(与长度L一致的矩形窗)
Px = Px / max(Px+realmin); % 归一到 1（便于阈值）
% (a)  -12 dB 平顶宽度
peak = 1; thr = db2mag(-12);
mask12 = (Px >= thr);
B12 = span_width(f, mask12);   % Hz
% (b)  95% 能量带宽（围绕谱峰从中心向外累计）
B95 = bandwidth95(Px, f);      % Hz
% 融合：取 max（偏保守，不漏边）
B_occ = max(B12, B95);
fprintf('估计占用带宽 B_occ ≈ %.3f MHz  [ -12dB: %.3f, 95%%能量: %.3f ]\n', B_occ/1e6, B12/1e6, B95/1e6);

% 展示 PSD（可选）
figure('Color','w','Name','PSD (norm)');
plot(f/1e6, 10*log10(Px+realmin),'LineWidth',1.0); grid on;
xlabel('Frequency (MHz)'); ylabel('PSD (dB, norm)');
title(sprintf('Normalized PSD, B_{occ}≈%.3f MHz', B_occ/1e6));

%% ===== 4) 由 B_occ 生成 Fs_guess 候选 =====
Fs_cands = B_occ ./ eta_grid;     % Hz
Fs_cands = unique(round(Fs_cands));  % 取整 Hz（也可 round 到 1e3/1e6）
fprintf('\nFs_guess 候选（由 η∈{%s} 推得）:\n', num2str(eta_grid,'%0.2f '));
disp(table(Fs_cands(:), 'VariableNames', {'Fs_guess_Hz'}));

%% ===== 5) 用 CP-相关峰做“轻验证”，给每个 Fs_guess 打分 =====
% 先一次 FFT-自相关，后续复用
[Rabs_full, taus_full] = r0_curve_abs_fft(y, 1, min(floor(L/2), round(Fr*0.1))); % τ 到 0.1s 上限或 L/2
% 为每个 Fs_guess + 每个 N，在 τ≈ N*(Fr/Fs_guess) 的 ±p 窗内看 |R0| 峰与动态范围
score_tab = [];  % 收集 [Fs_guess, best_dyn_dB, best_N, best_tau]
for j = 1:numel(Fs_cands)
    Fs_g = Fs_cands(j);
    best_dyn = -inf; best_N = NaN; best_tau = NaN; best_peak = -inf;
    for k = 1:numel(S_N)
        N = S_N(k);
        tau_c = N * (Fr / Fs_g);           % 期望滞后（样点）
        if tau_c < 2 || tau_c > taus_full(end), continue; end
        lo = max(1, floor(tau_c*(1-p)));
        hi = min(taus_full(end), ceil(tau_c*(1+p)));
        mask = (taus_full>=lo & taus_full<=hi);
        if ~any(mask), continue; end
        seg = Rabs_full(mask);
        pk = max(seg);
        dyn = 20*log10( max(seg)/max(min(seg), realmin) );
        % 评分：以 dyn 为主，峰值次之
        sc = dyn + 0.5*mag2db(pk+realmin);
        if sc > (best_dyn + 0.5*mag2db(best_peak+realmin))
            best_dyn = dyn; best_peak = pk; best_N = N;
            % 对应 τ 取峰位置
            [~,ix] = max(seg);
            taus_win = taus_full(mask);
            best_tau = taus_win(ix);
        end
    end
    score_tab = [score_tab; Fs_g, best_dyn, best_N, best_tau]; %#ok<AGROW>
end
score_tbl = array2table(score_tab, ...
    'VariableNames', {'Fs_guess_Hz','dyn_dB','N_best','tau_best_samp'});

% 排序：dyn_dB 降序；若接近则峰值更大者已在上面融进去
score_tbl = sortrows(score_tbl, 'dyn_dB', 'descend');

%% ===== 6) 输出推荐 Fs_guess =====
fprintf('\n—— 候选 Fs_guess 按“动态范围”打分排序（越大越可信）——\n');
disp(score_tbl);

rec = score_tbl(1,:);
fprintf('推荐 Fs_guess = %.6f MHz  [N_best=%d, τ_best≈%.0f samp, dyn≈%.1f dB]\n', ...
    rec.Fs_guess_Hz/1e6, rec.N_best, rec.tau_best_samp, rec.dyn_dB);

% 同时给出反推 η_hat 供自检
eta_hat = B_occ / rec.Fs_guess_Hz;
fprintf('回验占用率 η_hat = %.3f （常见范围 0.80~0.95）\n', eta_hat);

% 小提示
if rec.dyn_dB < nu_soft
    fprintf('提示：最佳 dyn 仅 %.1f dB，可能 SNR 较低/训练段干扰；可增大 p 或更换稳态片段重试。\n', rec.dyn_dB);
end

end % main

%% ================== 辅助函数 ==================
function raw = read_int16(filename)
    fid = fopen(filename,'rb','ieee-le'); assert(fid>0,'无法打开文件: %s',filename);
    c = onCleanup(@() fclose(fid));
    raw = fread(fid, inf, 'int16=>int16');
end

function w = span_width(f, mask)
    % 由布尔掩码得到连续覆盖频带宽度（Hz）
    if ~any(mask), w = 0; return; end
    idx = find(mask);
    w = abs(f(idx(end)) - f(idx(1)));
end

function B = bandwidth95(Px_norm, f)
    % 围绕谱峰，从峰向两侧扩展，直到累计功率达到 95% 的频宽
    [~,ix] = max(Px_norm);
    L = numel(Px_norm);
    left = ix; right = ix; s = Px_norm(ix);
    total = sum(Px_norm);
    while s < 0.95*total && (left>1 || right<L)
        % 按“哪侧下一个点更大就先取哪侧”的贪心扩展
        candL = left-1; candR = right+1;
        valL = (candL>=1) * Px_norm(max(candL,1));
        valR = (candR<=L) * Px_norm(min(candR,L));
        if valL >= valR
            if candL>=1, left=candL; s=s+Px_norm(left); else, right=candR; s=s+Px_norm(right); end
        else
            if candR<=L, right=candR; s=s+Px_norm(right); else, left=candL; s=s+Px_norm(left); end
        end
    end
    B = abs(f(right) - f(left));
end

function [Rabs, taus] = r0_curve_abs_fft(y, tmin, tmax)
    % 线性自相关 |R0(tau)|，一次 FFT 算所有 τ，再按 (L-τ) 归一
    y = y(:);
    L = numel(y);
    Nfft = 2^nextpow2(L + tmax);
    Y = fft(y, Nfft);
    R = ifft(abs(Y).^2);     % 线性相关的前 (tmax) 有效
    taus = tmin:tmax;
    Rsel = R(taus);
    m = (L - taus).';
    Rabs = abs(Rsel) ./ max(m,1);
    Rabs(~isfinite(Rabs)) = 0;
end

function x = db2mag(db), x = 10.^(db/20); end
function y = mag2db(x), y = 20*log10(abs(x)); end
