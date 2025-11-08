Fr = 409.6e6;     
% Fr = 240e6;% 把它换成你采数时的采样率
% opts.format     = 'int16_iq'; 
% opts.format     = 'float32_iq';% ✅ 正确类型
% opts.scale      = 1/32768;        % int16 归一化
opts.Detrend    = true;   
opts.Verbose    = true;    % 显示各候选过程
opts.NgTopK     = 10;      % 显示前 K 个 Ng 候选

% 设置多个 Fs_guess 候选值（单位：Hz）
opts.Fs_guess   = [320e6, 330e6, 340e6, 350e6, 334e6];  % 可根据需要调整候选值

% ===== CP 成对峰验证参数（新增）=====
opts.EnableCPCheck  = true;   % 启用 CP 检查（排除非 OFDM 周期）
opts.CPRatioMin     = 0.15;   % CP 副峰至少要达到主峰的 15%
opts.CPTolerancePct = 0.15;   % CP 位置容差 ±15%
opts.CPPenaltyDB    = 15;     % 没有 CP 副峰时惩罚 15 dB

R  = starlink_blind_id('test3.dat', Fr, opts);

% 查看所有候选结果
fprintf('\n========== 详细结果分析 ==========\n');
for i = 1:numel(R.all_Fs_candidates)
    res = R.all_Fs_candidates{i};
    fprintf('\n候选 #%d:\n', i);
    fprintf('  Fs_guess  = %.3f MHz\n', res.Fs_guess/1e6);
    fprintf('  N_hat     = %d\n', res.N_hat);
    fprintf('  Nr_hat    = %d\n', res.Nr_hat);
    fprintf('  Fs_hat    = %.6f MHz\n', res.Fs_hat/1e6);
    fprintf('  动态范围  = %.2f dB\n', res.diagN.dyn_dB_best);
    if isfield(res.diagN, 'cp_ratio_best')
        fprintf('  CP 副峰比 = %.1f%%\n', res.diagN.cp_ratio_best*100);
    end
end

fprintf('\n最终选择: 候选 #%d (Fs_guess=%.3f MHz)\n', ...
    R.best_candidate_idx, opts.Fs_guess(R.best_candidate_idx)/1e6);

% 可视化 R0 曲线（帮助诊断）
if isfield(R.diagN, 'taus_extended') && isfield(R.diagN, 'R0mag_extended')
    figure('Name', 'R0 自相关曲线 - CP 成对峰分析');
    plot(R.diagN.taus_extended, R.diagN.R0mag_extended, 'b-', 'LineWidth', 1.5);
    hold on;
    
    % 标注主峰
    plot(R.Nr_hat, interp1(R.diagN.taus_extended, R.diagN.R0mag_extended, R.Nr_hat), ...
        'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('主峰 @ τ=%d', R.Nr_hat));
    
    % 标注预期 CP 位置
    Ng_expected = round(R.Nr_hat * R.Ng_hat / (R.N_hat + R.Ng_hat));
    if Ng_expected > 0 && Ng_expected <= max(R.diagN.taus_extended)
        plot(Ng_expected, interp1(R.diagN.taus_extended, R.diagN.R0mag_extended, Ng_expected), ...
            'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('CP峰预期 @ τ=%d', Ng_expected));
    end
    
    xlabel('延迟 τ (样本)');
    ylabel('|R^0(τ)|');
    title(sprintf('自相关曲线 (N=%d, Ng=%d, Nr=%d)', R.N_hat, R.Ng_hat, R.Nr_hat));
    legend('Location', 'best');
    grid on;
    hold off;
end
