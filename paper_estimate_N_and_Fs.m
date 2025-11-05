function [N_hat, Nr_hat, Fs_hat, diagN] = paper_estimate_N_and_Fs(y, Fr, S_N, Fs_guess, p, nu_dB, M, Verbose)
% 构造接收机域区间 Sr（式(17)-(19)），扫 |R0(τ)|（式(11)近似），做 ν(dB) 验证（式(20)），得 N、Nr、Fs（式(22)）
y = y(:); L = numel(y);
if nargin<7 || isempty(M), M = min(L-1, max(round(0.015*Fr), round(0.5*L))); end
if nargin<8 || isempty(Verbose), Verbose = false; end
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

% 分桶找峰 + 每候选的 dyn(dB)
peak_vals = -inf(numel(S_N),1);
tau_arg   = zeros(numel(S_N),1);
dyn_vec   = -inf(numel(S_N),1);
for k=1:numel(S_N)
    tset=cells{k}; if isempty(tset), continue; end
    mask = ismember(taus,tset);
    seg  = R0mag(mask);
    [mx,ix]=max(seg);
    peak_vals(k)=mx; tau_arg(k)=tset(ix);
    dyn_vec(k) = 20*log10(max(seg)/max(min(seg),eps));
    if Verbose
        passflag = dyn_vec(k) >= nu_dB;
        fprintf('[N/Fs] N=%d: tau=[%d,%d], peak=%.4g@%d, dyn=%.1f dB%s\n', ...
            S_N(k), tset(1), tset(end), mx, tau_arg(k), dyn_vec(k), ternary_str(passflag,' PASS',''));
    end
end

% 验证门限：删掉不合格者（式(20)）
peak_vals_f = peak_vals;
while true
    [~,i_best]=max(peak_vals_f);
    if ~isfinite(peak_vals_f(i_best)), error('没有通过门限的候选；请增大 p 或 M。'); end
    if dyn_vec(i_best) >= nu_dB, dyn_dB_best = dyn_vec(i_best); break; else, peak_vals_f(i_best)=-inf; end
end

N_hat  = S_N(i_best);
Nr_hat = tau_arg(i_best);
Fs_hat = round( N_hat * (Fr / max(Nr_hat,1)) ); % 式(22)

if Verbose
    fprintf('[N/Fs] Chosen: N=%d, Nr=%d, Fs_hat=%.6f MHz, dyn=%.1f dB\n', N_hat, Nr_hat, Fs_hat/1e6, dyn_dB_best);
end

diagN = struct('taus',taus,'R0mag',R0mag,'S',S_N,'peak_vals',peak_vals,'tau_arg',tau_arg,'dyn_dB',dyn_vec,'dyn_dB_best',dyn_dB_best);
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
