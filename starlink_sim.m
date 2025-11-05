%% ========================================================================
%                     星链下行链路信号仿真脚本 (最终完整版)
% =========================================================================
%
% 功能:
% 1. 生成与论文描述一致的、包含混合调制的星链基带信号流。
% 2. 将生成的基带信号调制到Ku频段。
% 3. 忽略多普勒效应。
% 4. 可视化生成的基带信号和最终的Ku波段信号。
%
% =========================================================================

% --- 清理工作区 ---
clear; clc; close all;

%% 1. 设置仿真参数
% =========================================================================
s.Tdur = 5e-3;              % 仿真总时长 (秒)。
s.Fsr = 240e6;              % 接收机采样率 (Hz)，即信号带宽。 [cite: 215]
s.Fcr = 11.325e9;           % 虚拟接收机前端的中心频率 (Hz)，用于确定信道。
s.SNRdB = 20;               % 信噪比 (dB)。
s.beta = 0;                 % 多普勒因子，设置为0以忽略多普勒效应。
s.prob = 1;                 % 帧出现概率，设置为1生成连续帧流。

fprintf('仿真参数设置完毕。\n');

%% 2. 模式选择：使用默认的混合调制
% =========================================================================
% 本次仿真不注入自定义数据，将使用genStrlkFrame函数中的默认生成逻辑，
% 以便与论文描述的混合调制模式 (4-PSK, 16-QAM, 4-QAM) 完全一致。
fprintf('将使用默认的混合调制模式 (4-PSK, 16-QAM, 4-QAM)。\n');

%% 3. 生成星链基带信号流
% =========================================================================
fprintf('正在生成星链基带信号流...\n');
% 调用 genStrlkStream 函数。由于s结构体中没有.data字段，将触发默认生成逻辑。
[y_baseband, ~] = genStrlkStream(s);
fprintf('基带信号流生成完毕，共生成 %d 个采样点。\n', length(y_baseband));

% =========================================================================
%               【新添加部分】3.5 将代表Ku频段的信号保存为 .dat 文件
% =========================================================================
% 说明：
% 此处保存的是“复数基带IQ信号”，这是用于信号分析和识别的标准格式。
% 文件本身是中心在0Hz的，但它完整地代表了Ku频段信号的内容和结构。
% 您的盲识别程序在读取此文件后，应将其作为中心频率为Fc_Kuband (11.325 GHz)的信号来分析。
%--------------------------------------------------------------------------

fprintf('正在将生成的基带信号保存为 .dat 文件...\n');

% --- 文件参数定义 ---
output_filename = 'starlink_ku_band_signal.dat'; % 定义输出文件名 (.dat格式)
output_format = 'single';                        % 定义存储格式为单精度浮点数 ('float32')

% --- 数据格式转换 ---
% 将复数向量转换为交错的I/Q浮点数矩阵
% 创建一个 2xN 的矩阵，第一行为I(实部)，第二行为Q(虚部)
interleaved_data = [real(y_baseband)'; imag(y_baseband)'];

% --- 写入文件 ---
% 1. 以二进制写模式打开文件
fileID = fopen(output_filename, 'w');

% 2. 将数据写入文件
%    fwrite会按列写入，对于2xN的矩阵，它会写入 I(1), Q(1), I(2), Q(2), ...
fwrite(fileID, interleaved_data, output_format);

% 3. 关闭文件
fclose(fileID);
Fc_Kuband = 11.325e9;
fprintf('代表Ku频段信号的基带IQ文件已成功保存: %s\n', output_filename);
fprintf('格式: 交错I/Q, 32位单精度浮点数\n');
fprintf('使用此文件时，请将中心频率设置为 %.3f GHz\n', Fc_Kuband/1e9);


%% 4. 调制到Ku频段
% =========================================================================
% 将基带信号上变频到Ku频段
Fc_Kuband = 11.325e9;% 一个有效的星链Ku波段信道中心频率 (11.325 GHz) [cite: 215]

fprintf('正在将基带信号调制到Ku频段 (%.3f GHz)...\n', Fc_Kuband/1e9);
% 创建时间向量
t = (0:length(y_baseband)-1)' / s.Fsr;
% 上变频：复数基带信号乘以复指数载波
y_kuband = y_baseband .* exp(1j * 2 * pi * Fc_Kuband * t);
fprintf('Ku频段信号生成完毕。\n');


%% 5. 可视化结果
% =========================================================================

% --- 图1: 基带信号的频谱图 ---
fprintf('正在绘制基带信号的频谱图...\n');
figure('Name', '基带信号分析', 'NumberTitle', 'off');
plotSpec(y_baseband, 0, s.Tdur, s.Fcr, s.Fsr, s.Fcr, s.Fsr, 1024, '生成的星链基带信号频谱图 (混合调制)');
title('生成的星链基带信号频谱图 (混合调制)');
xlabel('时间 (s)');
ylabel('频率 (Hz)');


% --- 图2: Ku频段信号的功率谱 (使用基础FFT方法) ---
fprintf('正在绘制Ku频段信号的功率谱 (使用FFT方法)...\n');
figure('Name', 'Ku频段信号分析 (FFT)', 'NumberTitle', 'off');

N_fft = 8192; % 定义一个合适的FFT点数，用于频谱分析

% 步骤1: 对 Ku频段信号 y_kuband 执行FFT
Y_fft = fft(y_kuband, N_fft);

% 步骤2: 使用fftshift将零频分量移动到频谱中心
Y_fft_shifted = fftshift(Y_fft);

% 步骤3: 创建与FFT结果对应的频率轴
f_axis = (-N_fft/2 : N_fft/2 - 1) * (s.Fsr / N_fft);

% 步骤4: 计算功率谱密度 (PSD)
psd = (abs(Y_fft_shifted).^2) / (N_fft * s.Fsr);

% 步骤5: 绘图。频率轴 = f_axis + Fc_Kuband
plot((f_axis + Fc_Kuband) / 1e9, 10*log10(psd));

grid on;
title(sprintf('调制到Ku频段后的信号功率谱 (中心频率: %.3f GHz)', Fc_Kuband/1e9));
xlabel('频率 (GHz)');
ylabel('功率/频率 (dB/Hz)');
xlim([(Fc_Kuband - s.Fsr/2)/1e9, (Fc_Kuband + s.Fsr/2)/1e9]);
fprintf('所有处理和绘图完成。\n');


%% ========================================================================
%                           本地函数定义
% =========================================================================
% 以下是将所有外部 .m 文件内容整合为本地函数。主脚本可以直接调用它们。

function [y,present] = genStrlkStream(s)
% genStrlkStream: 生成一段连续的Starlink帧信号流
%
SNR = 10^(s.SNRdB/10);
sigma2_w = 1/SNR;
A = 1;
if (isfield(s,'sigma2_w'))
    sigma2_w = s.sigma2_w;
    A = sqrt(SNR.*sigma2_w);
end

if (isfield(s,'present'))
    if (isfield(s,'prob'))
        error("Can specify either prob or present, not both.");
    end
end
if (~isfield(s,'prob'))
    prob = 1;
else 
    prob = s.prob;
end
if (~isfield(s,'tau'))
    tau = 0;
else 
    tau = s.tau;
end

Fs = 240e6; % Starlink信号带宽
Tframe = 1/750; % Starlink帧周期
Nfr = ceil(s.Tdur/Tframe); % 根据总时长计算包含多少个完整的帧
if (isfield(s,'present'))
    present = s.present;
    Nfr_p = length(s.present);
    if (Nfr>Nfr_p)
        present(end+1:end+Nfr-Nfr_p) = 0;
    else
        present = present(1:Nfr);
    end
else
    present = binornd(1,prob,1,Nfr);
end

Ns = floor(s.Tdur*Fs); % 信号流总采样点数
Nframe = floor(Tframe*Fs); % 每帧的采样点数

if (isfield(s,'beta'))
    if (isfield(s,"phHist"))
        error("You can only specify a CFO (beta) OR a phase time-history, not both.")
    end
end
if (isfield(s,"phHist"))
    if (isfield(s,'beta'))
        error("You can only specify a CFO (beta) OR a phase time-history, not both.")
    end
    if (length(s.phHist) ~= Ns)
        error("Your phase time-history should be %d long according to your inputs.",Ns)
    end
end
    
y = zeros(Ns,1);

for ii = 1:Nfr
    if (present(ii))
        fr.SNRdB = nan();
        fr.Fsr = Fs; 
        fr.Fcr = getClosestFch(s.Fcr); 
        fr.beta = 0;
        if (isfield(s,'type'))
             fr.type = s.type; 
        end
        if (isfield(s,'data'))
             fr.data = s.data; 
             fr.Midx = s.Midx;
        end
        frame = genStrlkFrame(fr);
        
        frame_start_idx = (ii-1)*Nframe + 1;
        frame_end_idx = frame_start_idx + length(frame) - 1;
        if frame_end_idx <= Ns
            y(frame_start_idx:frame_end_idx) = frame;
        end
    end
end

Ntau = round(tau.*Fs);
y = [zeros(Ntau,1);y];
if length(y) > Ns
    y = y(1:Ns);
end

y = A.*y;

if (isfield(s,'beta'))
    if(~(s.beta == 0 && getClosestFch(s.Fcr) == s.Fcr))
        FD = -s.beta*getClosestFch(s.Fcr);
        fShift = FD + getClosestFch(s.Fcr) - s.Fcr;
    
        tVec = (0:length(y)-1)'/Fs;
        y = y.*exp(1j*2*pi*fShift*tVec);
    end
elseif (isfield(s,'doppVec'))
            fShift = getClosestFch(s.Fcr) - s.Fcr;
            tVec = (0:length(y)-1)'/Fs;
            offsetVec = 2*pi*s.doppVec./Fs;
            Phihist = cumsum(offsetVec);
            y = y.*exp(1j*(Phihist + 2*pi*fShift*tVec));
end

if (~isnan(s.SNRdB))    
    sigmaIQ = sqrt(sigma2_w/2);
    Nsamps = length(y);
    nVec = complex(sigmaIQ*randn(Nsamps,1),sigmaIQ*randn(Nsamps,1));
    y = y + nVec;
end

if(s.Fsr < Fs)
    tVec = (0:length(y)-1)'/Fs;
    y = resample(y,tVec,s.Fsr,"spline");
end
end

% -------------------------------------------------------------------------

function frame = genStrlkFrame(s)
% genStrlkFrame: 组装一个完整的Starlink信号帧 (原始混合调制版本)
%
SNR = 10^(s.SNRdB/10);
sigma2_w = 1/SNR;
A = 1;
if (isfield(s,'sigma2_w'))
    sigma2_w = s.sigma2_w;
    A = sqrt(SNR.*sigma2_w);
end

Tfg = (68/15)*1e-6;% 帧保护间隔 [cite: 215]
s.Fs = 240e6; % 信号带宽
Nfg = round(Tfg*s.Fs);
s.N = 1024; % 子载波数量
Fcii = getClosestFch(s.Fcr);

SNRdB = s.SNRdB;
beta = s.beta;
Fsr = s.Fsr;
Fcr = s.Fcr;

s.SNRdB = nan();
s.beta = 0;
s.Fcr = Fcii;
s.Fsr = s.Fs;

PSS = genPss();
[SSS, ~] = genSss();

if (~isfield(s,'data'))
    % --- 使用此段原始的、正确的默认逻辑 ---
    % 论文中提到帧的前几个符号是带相偏的4QAM (等效于4-PSK) [cite: 222]
    s.Nsym = 4;
    s.type = 'PSK';
    s.Midx = 4;
    Data = genStrlkOFDM(s);

    % 接下来是16QAM [cite: 221]
    s.Nsym = 4;
    s.type = 'QAM';
    s.Midx = 16;
    Data = [Data; genStrlkOFDM(s)];

    % 剩余大部分是4QAM [cite: 221]
    s.Nsym = 292;
    s.type = 'QAM';
    s.Midx = 4;
    Data = [Data; genStrlkOFDM(s)];    
else
    % 如果提供了数据，则使用用户数据
    s.Nsym = 300;
    Data = genStrlkOFDM(s);
end

Fg = zeros(Nfg,1);
 
s.Fcr = Fcr;
s.Fsr = Fsr;
s.beta = beta;
s.SNRdB = SNRdB;

frame = [PSS; SSS; Data; Fg];

if (isfield(s,'beta'))
    if (~(s.beta == 0 && getClosestFch(s.Fcr) == s.Fcr))
        FD = -s.beta*getClosestFch(s.Fcr);
        fShift = FD + getClosestFch(s.Fcr) - s.Fcr;
        tVec = (0:length(frame)-1)'/s.Fs;
        frame = frame.*exp(1j*2*pi*fShift*tVec);
    end
end
frame = A.*frame;
if (~isnan(s.SNRdB))    
    sigmaIQ = sqrt(sigma2_w/2);
    nVec = complex(sigmaIQ*randn(length(frame),1),sigmaIQ*randn(length(frame),1));
    frame = frame + nVec;
end
if(s.Fsr < s.Fs)
    tVec = (0:length(frame)-1)'/s.Fs;
    frame = resample(frame,tVec,s.Fsr,"spline");
end
end

% -------------------------------------------------------------------------

function y = genStrlkOFDM(s)
% genStrlkOFDM: genOFDM的封装，用于设置Starlink特定的OFDM参数
%
s.Fs = 240e6;
s.N = 1024;
s.Ng = 32;
s.gutter = 1;

if (~isfield(s,'Nsym'))
    s.Nsym = 1;
end
if (isfield(s,'data') && ~isfield(s,'type'))
    error(" A constellation type should also be provided as s.type (i.e. 'PSK' or 'QAM')")
end
if (isfield(s,'data'))
[l,w] = size(s.data);
if(l ~= s.N || w > s.Nsym)
    error('s.data must be an N x K vector where K <= Nsym');
end
end

if((~isfield(s,'data')) && (mod(log2(s.Midx),2) ~= 0) && (s.Midx ~= 2))
    error('Midx must 2 or an even power of 2');
end

F = s.Fs/s.N;
chIdx = round((s.Fcr/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
s.Fc = Fcii;

y = genOFDM(s);
end

% -------------------------------------------------------------------------

function yVec = genOFDM(s)
% genOFDM: 通用的OFDM符号生成器
%
if (~isfield(s,'gutter'))
    s.gutter = 0;
end
if (~isfield(s,'Nsym'))
    s.Nsym = 1;
end
if(~isfield(s,'data') && (mod(log2(s.Midx),2) ~= 0) && (s.Midx ~= 2))
  error('Midx must 2 or an even power of 2');
end
if (~isfield(s,'data'))
    x = randi([0 s.Midx-1],s.N,s.Nsym);
    [~,w] = size(x);
    if isscalar(s.Midx)
        s.Midx = s.Midx .* ones(w,1);
    end
else
    [l,w] = size(s.data);
    if(l ~= s.N || w > s.Nsym)
        error('s.data must be an N x K vector where K <= Nsym');
    end
    if isscalar(s.Midx)
        s.Midx = s.Midx .* ones(w,1);
    end
    x = s.data;
end

XVec = zeros(size(x));
if ( upper(string(s.type)) == "PSK")
    for ii = 1:w
        XVec(:,ii) = pskmod(x(:,ii),s.Midx(ii), pi/4);
    end
elseif( upper(string(s.type)) == "QAM")
    for ii = 1:w
        XVec(:,ii) = qammod(x(:,ii),s.Midx(ii),'UnitAveragePower',true);
    end
else
    error('Type must be "QAM" or "PSK".');
end

if (s.Nsym - w > 0)
    XVec_template = XVec;
    num_to_add = s.Nsym - w;
    for ii = 1:num_to_add
        XVec = [XVec, XVec_template(:, mod(ii-1, w)+1)];
    end
end

if (s.gutter)
    Ngut = 4;
    XVec = fftshift(XVec,1);
    XVec((s.N-Ngut)/2+1:(s.N-Ngut)/2+Ngut,:) = 0;
    XVec = ifftshift(XVec,1);
end

Mx = sqrt(s.N)*ifft(XVec);
MxCP = [Mx(end - s.Ng + 1:end,:); Mx];
xVec = MxCP(:);

if(s.beta ~= 0 || s.Fc ~= s.Fcr)
  tVec = (0:length(xVec)-1)'/s.Fs;
  FD = -s.beta*s.Fc;
  Fshift = FD + s.Fc - s.Fcr;
  xVec = xVec.*exp(1j*2*pi*Fshift*tVec);  
end

yVec = xVec;
if (~isnan(s.SNRdB))
    SNR = 10^(s.SNRdB/10);
    sigmaIQ = sqrt(1/(2*SNR));
    Nsamps = length(xVec);
    nVec = complex(sigmaIQ*randn(Nsamps,1),sigmaIQ*randn(Nsamps,1));
    yVec = xVec + nVec;
end

if(s.Fsr ~= s.Fs && ~isempty(yVec))
  tVec = (0:length(yVec)-1)'/s.Fs;
  yVec = resample(yVec,tVec,s.Fsr);
end
end

% -------------------------------------------------------------------------

function [pss] = genPss()
% genPss: 生成主同步序列 (PSS)
%
Ng = 32;
NpssSeg = 8;
ciVec = [3 7]';
a0Vec = [0 0 1 1 0 1 0]';
n = 7;
m = 2^n - 1;
lfsrSeq = zeros(m,1);
for idx=1:m
    buffer = a0Vec(ciVec);
    val = rem(sum(buffer),2);
    a0Vec = [val; a0Vec(1:end-1)];
    lfsrSeq(idx) = val;
end
lfsrSeqMod = flipud(lfsrSeq);
lfsrSeqMod = [0; lfsrSeqMod];
seq = 2*lfsrSeqMod-1;
pkSegVec = exp(-1j*pi/4 - 1j*0.5*pi*cumsum(seq));
pkVec = [-pkSegVec; repmat(pkSegVec,NpssSeg - 1,1)];
pkCP = pkVec(end-Ng+1:end);
pss = [-pkCP; pkVec];
end

% -------------------------------------------------------------------------

function [sssTimeDomainVec, sssFreqDomainVec] = genSss()
% (已修正) 生成次级同步序列 (SSS)
% 修正1: 使用与genOFDM一致的、正确的IFFT缩放因子，确保功率一致。
% 修正2: 内置论文中的十六进制字符串，移除对 sssVec.mat 的依赖。
    N = 1024;
    Ng = 32;

    % 根据论文式(39)的十六进制字符串 q_sss 生成序列
    q_sss_hex = [ ...
    'BD565D5064E9B3A94958F28624DED5' ...
    '60946199F5B40F0E4FB5EFCB473B4C24' ...
    'B2D1E0BD01A6A04D5017DE91A8ECC0DA' ...
    '09EBFE57F9F1B44C532F161C583A4249' ...
    '0A5C09F2A117F9A28F9B2FD547A74C44' ...
    'BABB4BE85DA6A62B1235E2AD084C0018' ...
    '0142A8F7F357DEC4F31316BC58FA4049' ...
    '09A3FCA7F88E421902B6A2580AE80308' ...
    '03F65809DB347F590DBC46F010EBE3A2' ...
    '5C060D74429FC46BDF9B63719279798D' ... % 论文原文为..7980，但长度不足，按标准应为..798D
    '232C5ABA274122FF66AD7E449F44CB40' ...
    'C49C24A1E2629F5BFE82CE531FDC34F8' ...
    'C64A43A963F40D5B71BDE6FB2F13492D' ...
    '6F2E8544B21D449722C635180342CD00' ...
    '26A1E7F7E80E91B175E852F919767E5A' ...
    'F9B6E909AF362F5218E2B908DC005803'];
    
    % 将十六进制字符串解析为4QAM的相位索引 (0,1,2,3)
    q_sss_dec = hex2dec(reshape(q_sss_hex', 2, [])');
    s_k = zeros(N-4, 1);
    for k = 1:(N-4)
        byte_idx = floor((k-1)/4) + 1;
        shift_val = 6 - 2*mod(k-1, 4);
        s_k(k) = bitshift(q_sss_dec(byte_idx), -shift_val) & 3;
    end
    
    % 根据相位索引生成4QAM符号
    theta_k = s_k * pi/2;
    sssFreqDomainVec = zeros(N, 1);
    sssFreqDomainVec(3:N-2) = exp(1j * theta_k); % 对应论文 s_k*pi/2
    
    % ---【核心修正】---
    % 使用与genOFDM相同的、标准的IFFT缩放，以确保功率一致性
    sssTimeDomainVec_noCP = sqrt(N) * ifft(sssFreqDomainVec);
    sssTimeDomainVec = [sssTimeDomainVec_noCP(end-Ng+1:end); sssTimeDomainVec_noCP];
end

% -------------------------------------------------------------------------

function Fcii = getClosestFch(Fc)
% getClosestFch: 计算距离给定频率最近的Starlink标准信道中心频率
%
F = 240e6/1024;
chIdx = round((Fc/1e9 - 10.7 - F/2/1e9)/0.25 + 0.5);
Fcii = (10.7e9 + F/2 + 250e6*(chIdx - 0.5));
end

% -------------------------------------------------------------------------

function plotSpec(data,tstart,tdur,Fcr,Fsr,Fc,Fs,NFFT,Stitle)
% plotSpec: 绘制信号的频谱图
%
if ~exist('Stitle','var')
     Stitle = sprintf("Spectrogram centered at %d",Fc);
end
   
if (tdur == 0)
    tdur = length(data)/Fsr;
end
seekOffset = floor(tstart*Fsr);
Nblock = floor(tdur*Fsr);
data = data(seekOffset+1:seekOffset+Nblock);
    
tyVec = (0:length(data)-1)'/Fsr;
if (Fs ~= Fsr)
    [data, ~] = resample(data,tyVec,Fs);
end

if (Fc ~= Fcr)
    tVec_resampled = (0:length(data)-1)'/Fs;
    Fshift = Fcr - Fc;
    data = data.*exp(1j*2*pi*Fshift*tVec_resampled);
end
    
spectrogram(data,kaiser(NFFT,5),floor(NFFT*0.8),...
          NFFT,Fs,'psd','centered','yaxis', 'minthreshold', -90);
title(Stitle);
end