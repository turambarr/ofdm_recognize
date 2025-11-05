function main_detect_ofdm_float
    % --- 1. 用户必须配置的参数 ---
    FILE_PATH = 'starlink_ku_band_signal.dat'; % 替换为您的float32 IQ文件路径
    SAMPLE_RATE = 240e6; % 真实采样率 (Hz)
    SAMPLES_TO_READ = 4e6; % 要读取和分析的采样点总数 (复样本数)

    % 功率谱图的起止样点 (设置为0则使用全部数据)
    PSD_START_SAMPLE = 1; % 起始样点索引 (从1开始，0表示从头开始)
    PSD_END_SAMPLE = 1200000;   % 结束样点索引 (0表示到末尾)
    
    % --- 2. 读取 I/Q 数据 (float32交织) ---
    fprintf('正在读取 %d 个 I/Q 样本从 %s (float32)...\n', SAMPLES_TO_READ, FILE_PATH);
    try
        iq_data = read_iq_dat_float(FILE_PATH, SAMPLES_TO_READ);
        fprintf('数据读取完毕。\n');
    catch e
        fprintf('读取文件时出错: %s\n', e.message);
        return;
    end
    
    % --- 3. 绘制功率谱密度 (PSD) ---
    fprintf('正在计算并绘制 PSD...\n');
    
    % 确定PSD绘制的样点范围
    if PSD_START_SAMPLE == 0
        psd_start = 1;
    else
        psd_start = PSD_START_SAMPLE;
    end
    
    if PSD_END_SAMPLE == 0
        psd_end = length(iq_data);
    else
        psd_end = min(PSD_END_SAMPLE, length(iq_data));
    end
    
    % 验证范围
    if psd_start < 1 || psd_start > length(iq_data)
        error('PSD起始样点超出范围！');
    end
    if psd_end < psd_start
        error('PSD结束样点必须大于起始样点！');
    end
    
    % 提取用于PSD的数据段
    iq_data_psd = iq_data(psd_start:psd_end);
    fprintf('PSD使用样点范围: %d 到 %d (共 %d 个样本)\n', ...
            psd_start, psd_end, length(iq_data_psd));
    
    figure(1);
    pwelch(iq_data_psd, hann(8192), 4096, 8192, SAMPLE_RATE, 'centered');
    title(sprintf('功率谱密度 (PSD) - 样本范围 [%d:%d]', psd_start, psd_end));
    xlabel('频率 (Hz)');
    ylabel('功率/频率 (dB/Hz)');
    
    fprintf('PSD 分析完成。\n');
end

% -------------------------------------------------------------------------
% 辅助函数: 读取 float32 I/Q .dat 文件 (I,Q交错)
% -------------------------------------------------------------------------
function iq_samples = read_iq_dat_float(file_path, num_samples)
    % 打开文件
    fid = fopen(file_path, 'r');
    if fid == -1
        error('无法打开文件: %s', file_path);
    end
    
    % 读取 single 数据，I 和 Q 是交错的，所以要读取 2 * num_samples 个值
    data_raw = fread(fid, num_samples * 2, 'single');
    
    % 关闭文件
    fclose(fid);
    
    if length(data_raw) < num_samples * 2
        warning('读取的样本数少于请求的样本数。');
    end
    
    % 分离 I 和 Q
    data_I = data_raw(1:2:end);
    data_Q = data_raw(2:2:end);
    
    % 合成复数信号（转 double 便于后续DSP与绘图）
    iq_samples = complex(double(data_I), double(data_Q));
    
    % 确保为列向量
    iq_samples = iq_samples(:);
end
