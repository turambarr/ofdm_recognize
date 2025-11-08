%% OFDM 盲识别 - CP 成对峰验证对比测试
% 对比启用/禁用 CP 检查的效果

Fr = 409.6e6;
opts.Detrend = true;
opts.Verbose = true;
opts.NgTopK = 10;

% Fs_guess 候选（根据实际调整）
opts.Fs_guess = [320e6, 330e6, 340e6, 350e6, 334e6];

%% 场景1：禁用 CP 检查（原始算法）
fprintf('\n');
fprintf('====================================================\n');
fprintf('   场景1：禁用 CP 检查（仅基于动态范围 dyn）\n');
fprintf('====================================================\n\n');

opts.EnableCPCheck = false;
R_without_CP = starlink_blind_id('findN.dat', Fr, opts);

fprintf('\n--- 场景1 最终结果 ---\n');
fprintf('N = %d\n', R_without_CP.N_hat);
fprintf('Ng = %d\n', R_without_CP.Ng_hat);
fprintf('Fs = %.6f MHz\n', R_without_CP.Fs_hat/1e6);
fprintf('Nr = %d\n', R_without_CP.Nr_hat);

%% 场景2：启用 CP 检查（改进算法）
fprintf('\n\n');
fprintf('====================================================\n');
fprintf('   场景2：启用 CP 检查（排除虚假周期性）\n');
fprintf('====================================================\n\n');

opts.EnableCPCheck = true;
opts.CPRatioMin = 0.15;       % CP 副峰至少 15% 主峰强度
opts.CPTolerancePct = 0.15;   % CP 位置容差 ±15%
opts.CPPenaltyDB = 15;        % 无 CP 峰惩罚 15 dB

R_with_CP = starlink_blind_id('findN.dat', Fr, opts);

fprintf('\n--- 场景2 最终结果 ---\n');
fprintf('N = %d\n', R_with_CP.N_hat);
fprintf('Ng = %d\n', R_with_CP.Ng_hat);
fprintf('Fs = %.6f MHz\n', R_with_CP.Fs_hat/1e6);
fprintf('Nr = %d\n', R_with_CP.Nr_hat);
if isfield(R_with_CP.diagN, 'cp_ratio_best')
    fprintf('CP ratio = %.1f%%\n', R_with_CP.diagN.cp_ratio_best*100);
end

%% 对比结果
fprintf('\n\n');
fprintf('====================================================\n');
fprintf('                    结果对比\n');
fprintf('====================================================\n');
fprintf('参数          | 无 CP 检查      | 有 CP 检查      | 差异\n');
fprintf('-----------------------------------------------------------\n');
fprintf('N             | %-15d | %-15d | %s\n', R_without_CP.N_hat, R_with_CP.N_hat, ...
    ternary(R_without_CP.N_hat == R_with_CP.N_hat, '相同', '不同'));
fprintf('Ng            | %-15d | %-15d | %s\n', R_without_CP.Ng_hat, R_with_CP.Ng_hat, ...
    ternary(R_without_CP.Ng_hat == R_with_CP.Ng_hat, '相同', '不同'));
fprintf('Nr            | %-15d | %-15d | %s\n', R_without_CP.Nr_hat, R_with_CP.Nr_hat, ...
    ternary(R_without_CP.Nr_hat == R_with_CP.Nr_hat, '相同', '不同'));
fprintf('Fs (MHz)      | %-15.3f | %-15.3f | %s\n', R_without_CP.Fs_hat/1e6, R_with_CP.Fs_hat/1e6, ...
    ternary(abs(R_without_CP.Fs_hat - R_with_CP.Fs_hat) < 0.1e6, '相同', '不同'));
fprintf('Tsym (us)     | %-15.3f | %-15.3f | \n', R_without_CP.Tsym*1e6, R_with_CP.Tsym*1e6);
fprintf('F_sub (kHz)   | %-15.3f | %-15.3f | \n', R_without_CP.F/1e3, R_with_CP.F/1e3);

if R_without_CP.Nr_hat ~= R_with_CP.Nr_hat
    fprintf('\n✓ CP 检查改变了选择结果！\n');
    fprintf('  原算法选择了 τ=%d（可能是虚假峰）\n', R_without_CP.Nr_hat);
    fprintf('  CP 检查选择了 τ=%d（有 CP 副峰，更可靠）\n', R_with_CP.Nr_hat);
else
    fprintf('\n✓ 两种方法选择了相同的结果（信号质量好）\n');
end

%% 可视化对比
figure('Name', '对比：有/无 CP 检查', 'Position', [100 100 1200 500]);

% 左图：无 CP 检查
subplot(1,2,1);
if isfield(R_without_CP.diagN, 'taus_extended')
    plot(R_without_CP.diagN.taus_extended, R_without_CP.diagN.R0mag_extended, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(R_without_CP.Nr_hat, interp1(R_without_CP.diagN.taus_extended, ...
        R_without_CP.diagN.R0mag_extended, R_without_CP.Nr_hat), ...
        'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('选择 τ=%d', R_without_CP.Nr_hat));
    xlabel('延迟 τ (样本)');
    ylabel('|R^0(τ)|');
    title('场景1：无 CP 检查');
    legend('Location', 'best');
    grid on;
    hold off;
end

% 右图：有 CP 检查
subplot(1,2,2);
if isfield(R_with_CP.diagN, 'taus_extended')
    plot(R_with_CP.diagN.taus_extended, R_with_CP.diagN.R0mag_extended, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(R_with_CP.Nr_hat, interp1(R_with_CP.diagN.taus_extended, ...
        R_with_CP.diagN.R0mag_extended, R_with_CP.Nr_hat), ...
        'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('主峰 τ=%d', R_with_CP.Nr_hat));
    
    % 标注 CP 副峰位置
    Ng_exp = round(R_with_CP.Nr_hat * R_with_CP.Ng_hat / (R_with_CP.N_hat + R_with_CP.Ng_hat));
    if Ng_exp > 0 && Ng_exp <= max(R_with_CP.diagN.taus_extended)
        plot(Ng_exp, interp1(R_with_CP.diagN.taus_extended, ...
            R_with_CP.diagN.R0mag_extended, Ng_exp), ...
            'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('CP峰 τ=%d', Ng_exp));
    end
    
    xlabel('延迟 τ (样本)');
    ylabel('|R^0(τ)|');
    title('场景2：有 CP 检查');
    legend('Location', 'best');
    grid on;
    hold off;
end

fprintf('\n完成！\n');

function s = ternary(cond, t, f)
    if cond, s = t; else, s = f; end
end
