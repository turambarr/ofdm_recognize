function [y, meta] = read_iq_autodetect(fname, opts)
% 读取交错IQ并判别是 int16_iq 还是 float32_iq（小端默认）
% - 强制：opts.format='int16_iq' | 'float32_iq'
% - 自动：opts.format='auto' 或留空；用前缀小段双读评分择优
%
% 可选参数（opts）：
%   format   : 'auto'(默认) | 'int16_iq' | 'float32_iq'
%   endian   : 'ieee-le'(默认) | 'ieee-be'
%   order    : 'IQ'(默认) | 'QI'          % 若为 QI 将在读完后交换 I/Q
%   scale    : 仅 int16 有效的缩放（默认 1/32768）
%   maxProbe : 判别时用于评分的“复样点”上限（默认 5e5）
%   Verbose  : true/false 打印评分明细（默认 false）
%
% 输出：
%   y        : complex double 列向量
%   meta     : struct，含 detected_format/endian/order/scale 与评分

    % --- 读取参数 ---
    if nargin<2 || isempty(opts), opts = struct(); end
    fmt     = get_opt(opts,'format','auto');
    endian  = get_opt(opts,'endian','ieee-le');
    order   = upper(get_opt(opts,'order','IQ'));
    scale16 = get_opt(opts,'scale',1/32768);
    maxProbe= get_opt(opts,'maxProbe',5e5);     % 参与判别的“复样点”上限
    verbose = get_opt(opts,'Verbose',false);

    if order~="IQ" && order~="QI", error('opts.order 只能为 ''IQ'' 或 ''QI'''); end

    % --- 强制模式：直接整文件读取 ---
    if ~isempty(fmt) && ~strcmpi(fmt,'auto')
        y = read_full_as(fname, fmt, endian, scale16);
        if order=="QI", y = complex(imag(y), real(y)); end
        meta = struct('detected_format',sprintf('%s(forced)',lower(fmt)), ...
                      'endian',endian,'order',char(order),'scale',scale16);
        return;
    end

    % --- 自动判别：按 int16 / float32 双读一小段，打分择优 ---
    s16 = probe_stats(fname, 'int16',  endian, scale16, maxProbe);
    s32 = probe_stats(fname, 'single', endian, 1.0,     maxProbe);

    if verbose
        fprintf('[AUTO] int16: std=%.3g  circ=%.3g  corrRI=%.3g  ac1=%.3g  qstep=%.3g  bad=%d\n', ...
            s16.stdsum, s16.circ, s16.corrRI, s16.ac1, s16.qstep, s16.bad);
        fprintf('[AUTO] float: std=%.3g  circ=%.3g  corrRI=%.3g  ac1=%.3g            bad=%d\n', ...
            s32.stdsum, s32.circ, s32.corrRI, s32.ac1, s32.bad);
    end

    % 先用“坏样本”快速规则
    if s32.bad && ~s16.bad
        chosen = 'int16_iq';
    elseif s16.bad && ~s32.bad
        chosen = 'float32_iq';
    else
        % 强规则：Re/Im 相关显著区分（float32被当int16时 corrRI 常很大）
        if (s16.corrRI > 0.6 && s32.corrRI < 0.3)
            chosen = 'float32_iq';
        elseif (s32.corrRI > 0.6 && s16.corrRI < 0.3)
            chosen = 'int16_iq';
        else
            % 综合分数（stdsum 大、circ 小、corrRI 小、ac1 适中；int16 读法再加 qstep 奖惩）
            score16 = s16.stdsum * (1 - s16.corrRI) / (1 + 5*s16.circ) + 0.2*s16.ac1 ...
                      - 0.5*log10(max(s16.qstep,1e-12));    % qstep 越小越像 int16→float
            score32 = s32.stdsum * (1 - s32.corrRI) / (1 + 5*s32.circ) + 0.2*s32.ac1;
            chosen  = ternary(score16 >= score32, 'int16_iq', 'float32_iq');
        end
    end

    % --- 用判定格式“整文件读取” ---
    y = read_full_as(fname, chosen, endian, scale16);
    if order=="QI", y = complex(imag(y), real(y)); end

    % 输出元信息
    meta = struct('detected_format',[chosen '(auto)'], ...
                  'endian',endian,'order',char(order),'scale',scale16, ...
                  'score16',s16,'score32',s32);

end

% ====== 本文件内局部函数（不生成额外 .m） ======
function v = get_opt(S,name,defv)
    if isstruct(S) && isfield(S,name) && ~isempty(S.(name)), v=S.(name); else, v=defv; end
end
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end

function y = read_full_as(fname, fmt, endian, scale16)
    fid=fopen(fname,'rb',endian); assert(fid>0,'无法打开文件: %s',fname);
    switch lower(fmt)
        case {'int16_iq','int16'}
            r=fread(fid,inf,'int16'); fclose(fid);
            L=floor(numel(r)/2)*2; r=r(1:L);
            y = complex(double(r(1:2:end))*scale16, double(r(2:2:end))*scale16);
        case {'float32_iq','single','float32'}
            r=fread(fid,inf,'single'); fclose(fid);
            L=floor(numel(r)/2)*2; r=r(1:L);
            y = complex(double(r(1:2:end)), double(r(2:2:end)));
        otherwise
            fclose(fid); error('未知 format: %s（需 int16_iq/float32_iq）', fmt);
    end
    y = y(:);
end

function st = probe_stats(fname, dtype, endian, scl, maxCs)
% 读前缀 maxCs 个“复样点”等效的标量（2*maxCs 个）做评分
    st = struct('stdsum',0,'circ',Inf,'ac1',0,'corrRI',1,'qstep',Inf,'bad',true);
    fid=fopen(fname,'rb',endian); if fid<=0, return; end
    r=fread(fid, 2*min(maxCs,1e6), dtype); fclose(fid);    % 上限：≤50万复样点
    L=floor(numel(r)/2)*2; r=r(1:L); if L<4, return; end

    y = complex(double(r(1:2:end))*scl, double(r(2:2:end))*scl);
    if ~all(isfinite(real(y))) || ~all(isfinite(imag(y))), return; end

    % 去直流（不做RMS归一，保留尺度用于 stdsum）
    y = y - mean(y);
    re = real(y); im = imag(y);
    vr = mean(re.^2); vi = mean(im.^2);
    if vr<=0 || vi<=0, return; end

    st.stdsum = sqrt(vr) + sqrt(vi);                                         % 能量感知
    st.circ   = abs(mean(y.^2))/max(mean(abs(y).^2),eps);                     % 圆对称偏差
    st.ac1    = abs(mean( conj(y(1:end-1)).*y(2:end) ));                      % 一阶自相关模
    covRI     = mean( (re-mean(re)).*(im-mean(im)) );
    st.corrRI = abs( covRI / sqrt(vr*vi) );                                   % 实虚相关

    % 量化阶梯吻合度（仅在“按 int16 读法”时有意义；此处统一计算，供比较）
    % 原始 int16 值应 ≈ re/scale16、im/scale16 为整数；float32 错读时不会贴近整数
    q1 = abs(re/ max(scl,eps) - round(re/ max(scl,eps)));
    q2 = abs(im/ max(scl,eps) - round(im/ max(scl,eps)));
    st.qstep = median([q1; q2]);                                              % 越小越像 int16→float

    st.bad = false;
end
