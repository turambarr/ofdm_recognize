function main_detect_ofdm
% main_detect_ofdm.m
% 读取 int16 I/Q 交织数据，绘制 PSD，并自动估计占用带宽（-12 dB 平顶与 95%能量带宽）。
% 使用方式：直接运行 main_detect_ofdm（可修改下方“用户配置”区）

    %% ===== 1) 用户配置 =====
    FILE_PATH       = 'test3.dat';      % 数据文件路径（int16 I/Q 交织）
    SAMPLE_RATE     = 409.6e6;            % 采样率 Hz（用于频轴标定；如不确定可填写实际 Fr）
    SAMPLES_TO_READ = 4e6;              % 要读取的复样点数（I/Q合一后的样本数）
    PSD_START_SAMPLE= 1;                % PSD 起始样本索引 (从1开始; 0表示从开头)
    PSD_END_SAMPLE  = 0;                % PSD 结束样本索引 (0表示到末尾)
    USE_HAMMING     = true;             % PSD 用 Hamming 窗（泄漏更低）；false 则用矩形窗
    NFFT_MAX        = 2^18;             % PSD 的 nfft 上限（自动取不超过此值的 2^k）

    %% ===== 2) 读取 I/Q 数据 =====
    if ~isfile(FILE_PATH)
        error('找不到文件：%s', FILE_PATH);
    end
    fprintf('正在读取 %g 个 I/Q 复样点，自 %s ...\n', SAMPLES_TO_READ, FILE_PATH);
    iq_data = read_iq_dat(FILE_PATH, SAMPLES_TO_READ);
    L = numel(iq_data);
    fprintf('数据读取完毕，实际读取到 %.3f M 复样点。\n', L/1e6);

    if isempty(SAMPLE_RATE) || ~isfinite(SAMPLE_RATE) || SAMPLE_RATE<=0
        error('SAMPLE_RATE（采样率）不能为空或非法。');
    end

    %% ===== 3) 选择 PSD 片段 =====
    if PSD_START_SAMPLE == 0, PSD_START_SAMPLE = 1; end
    if PSD_END_SAMPLE   == 0, PSD_END_SAMPLE   = L; end
    if PSD_START_SAMPLE < 1 || PSD_START_SAMPLE > L
        error('PSD_START_SAMPLE 超出范围。');
    end
    if PSD_END_SAMPLE < PSD_START_SAMPLE
        error('PSD_END_SAMPLE 必须 >= PSD_START_SAMPLE。');
    end
    iq_psd = iq_data(PSD_START_SAMPLE:PSD_END_SAMPLE);
    Lseg   = numel(iq_psd);
    fprintf('PSD 使用样本范围: [%d : %d]，共 %d 样本。\n', PSD_START_SAMPLE, PSD_END_SAMPLE, Lseg);

    %% ===== 4) 计算 PSD（Welch）并估计占用带宽 =====
    x = iq_psd - mean(iq_psd);                 % 去 DC 偏置
    % 选择窗长：不超过 8192，且不超过片段长度
    winLen = min(Lseg, 8192);
    if USE_HAMMING
        win = hamming(winLen);
    else
        win = rectwin(winLen);                 % 矩形窗
    end
    nover = max(0, floor(winLen/2));
    nfft  = max(4096, 2^nextpow2(min(Lseg, NFFT_MAX)));

    [Px, f] = pwelch(x, win, nover, nfft, SAMPLE_RATE, 'centered'); % 单位：功率/Hz
    % 原先：Px_norm = Px / max(Px);
% 改为：相对“平顶”归一化（更稳健）
Px_sorted = sort(Px,'descend');
P_plateau = median(Px_sorted(round(0.01*numel(Px)):round(0.10*numel(Px)))); % 去掉最尖1%，取前10%的中位
Px_norm   = Px / max(P_plateau, realmin);

% 然后再用阈值法
B12 = bandwidth_thres(Px_norm, f, -12);  % 现在的 -12dB 是相对“平顶”                              % 归一到峰 = 1

    % (-12 dB 平顶) 与 (95% 能量) 两种带宽
    B12 = bandwidth_thres(Px_norm, f, -12);   % Hz
    B95 = bandwidth95(Px_norm, f);            % Hz
    Bocc = max(B12, B95);
% —— 估计“平顶”功率（去掉最尖1%，在前10%内取中位）——
Px_sorted = sort(Px,'descend');
i1 = max(1, round(0.01*numel(Px_sorted)));
i2 = max(i1+1, round(0.10*numel(Px_sorted)));
P_plateau = median(Px_sorted(i1:i2));

% —— 按“平顶”为0 dB 重新归一化，测 -3 / -6 dB 宽度 —— 
Px_rel = Px / max(P_plateau, realmin);
B3  = bandwidth_thres(Px_rel, f, -3);
B6  = bandwidth_thres(Px_rel, f, -6);

fprintf('平顶口径：B_-3dB = %.3f MHz,  B_-6dB = %.3f MHz\n', B3/1e6, B6/1e6);
B3_peak = bandwidth_thres(Px_norm, f, -3);
B6_peak = bandwidth_thres(Px_norm, f, -6);
fprintf('相对峰值：B_-3dB(peak) = %.3f MHz,  B_-6dB(peak) = %.3f MHz\n', B3_peak/1e6, B6_peak/1e6);

    %% ===== 5) 画图展示 =====
    figure('Color','w','Name','PSD (Welch, centered)');
    plot(f/1e6, 10*log10(Px + realmin), 'LineWidth',1.05); grid on;
    xlabel('Frequency (MHz)');
    ylabel('PSD (dB/Hz)');
    title(sprintf('Welch PSD [%d:%d], win=%s, nfft=%d', ...
        PSD_START_SAMPLE, PSD_END_SAMPLE, ternary(USE_HAMMING,'Hamming','Rect'), nfft));

    % 标注 -12 dB 覆盖区（可视化）
    thr = db2mag(-12);
    mask12 = (Px_norm >= thr);
    yl = ylim;
    hold on;
    if any(mask12)
        plot(f(mask12)/1e6, repmat(yl(1)+1, sum(mask12),1), '.', 'Color', [0.85 0.2 0.2]);
    end
% 画图后，先取坐标轴范围
xLim = get(gca,'XLim');
yLim = get(gca,'YLim');

% 计算放文字的位置：x 轴左侧偏 2%，y 轴顶部往下 3 dB
xpos = xLim(1) + 0.02 * diff(xLim);
ypos = yLim(2) - 3;

text(xpos, ypos, ...
    sprintf('B_{12dB}=%.3f MHz,  B_{95%%}=%.3f MHz,  B_{occ}=%.3f MHz', ...
            B12/1e6, B95/1e6, Bocc/1e6), ...
    'FontWeight','bold');

    fprintf('占用带宽估计：B_12dB = %.3f MHz,  B_95%% = %.3f MHz,  取 B_{occ} = %.3f MHz\n', ...
        B12/1e6, B95/1e6, Bocc/1e6);

    fprintf('完成。提示：采样率 SAMPLE_RATE 仅用于设置频轴单位；实际“占用带宽”由 PSD 估计得到。\n');
end

%% ===================== 子函数们 =====================

function iq_samples = read_iq_dat(file_path, num_samples)
% 读取 int16 I/Q 交织数据；返回复基带列向量
    fid = fopen(file_path, 'rb', 'ieee-le');
    assert(fid>0, '无法打开文件: %s', file_path);
    c = onCleanup(@() fclose(fid));
    if isinf(num_samples) || num_samples<=0
        raw = fread(fid, inf, 'int16=>int16');
    else
        raw = fread(fid, 2*round(num_samples), 'int16=>int16');
    end
    if mod(numel(raw),2)~=0
        warning('原始样本数为奇数，丢弃最后1个样本以对齐 I/Q。');
        raw = raw(1:end-1);
    end
    I = double(raw(1:2:end));
    Q = double(raw(2:2:end));
    iq_samples = complex(I, Q);
    iq_samples = iq_samples(:);
end

function B = bandwidth_thres(Px_norm, f, dB)
% 在“相对幅度阈值（如 -12 dB）”处测量频宽
    thr = db2mag(dB);
    mask = (Px_norm >= thr);
    B = span_width(f, mask);
end

function B = bandwidth95(Pn, f)
% 从谱峰向两侧扩展，覆盖到 95% 累积能量的频宽
    [~,ix] = max(Pn);
    L = numel(Pn);
    left = ix; right = ix; s = Pn(ix); total = sum(Pn);
    while s < 0.95*total && (left>1 || right<L)
        candL = left-1;  candR = right+1;
        vL = (candL>=1) * Pn(max(candL,1));
        vR = (candR<=L) * Pn(min(candR,L));
        if vL >= vR
            if candL>=1, left=candL;  s=s+Pn(left);  else, right=candR; s=s+Pn(right); end
        else
            if candR<=L, right=candR; s=s+Pn(right); else, left=candL; s=s+Pn(left); end
        end
    end
    B = abs(f(right) - f(left));
end

function w = span_width(f, mask)
% 从布尔掩码得到连续覆盖的频带宽度（Hz）
    if ~any(mask), w = 0; return; end
    idx = find(mask);
    w = abs(f(idx(end)) - f(idx(1)));
end

function y = ternary(cond, a, b)
% 三元选择（字符串友好）
    if cond, y = a; else, y = b; end
end

function x = db2mag(db)
% dB → 幅度
    x = 10.^(db/20);
end
