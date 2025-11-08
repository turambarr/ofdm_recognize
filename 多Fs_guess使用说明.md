# 多 Fs_guess 候选值功能使用说明

## 功能概述

修改后的 `starlink_blind_id.m` 现在支持同时测试多个 `Fs_guess` 候选值，并自动选择最佳结果。

## 使用方法

### 1. 基本用法

```matlab
Fr = 409.6e6;  % 采样率

opts = struct();
opts.Detrend = true;
opts.Verbose = true;

% 设置多个 Fs_guess 候选值（数组形式）
opts.Fs_guess = [200e6, 240e6, 250e6, 300e6];

% 运行识别
R = starlink_blind_id('middletest.dat', Fr, opts);
```

### 2. 单个候选值（向后兼容）

如果只提供一个值，功能与之前完全相同：

```matlab
opts.Fs_guess = 240e6;  % 单个值
R = starlink_blind_id('middletest.dat', Fr, opts);
```

## 输出结果

### 返回结构 R 的新增字段

- `R.all_Fs_candidates`: cell 数组，包含所有候选的详细结果
  - 每个元素是一个结构体，包含：
    - `Fs_guess`: 该候选的 Fs_guess 值
    - `N_hat`: 估计的子载波数
    - `Nr_hat`: 估计的接收机域样本数
    - `Fs_hat`: 估计的信道带宽
    - `diagN`: 诊断信息（包括动态范围等）

- `R.best_candidate_idx`: 最佳候选的索引

### 选择最佳结果的标准

系统自动基于**动态范围 (dyn_dB)** 选择最佳候选。动态范围越高，表明该候选的估计质量越好。

## 示例输出

运行时会显示类似以下的输出：

```
========== 开始估计 N 和 Fs（共 4 个 Fs_guess 候选）==========

---------- Fs_guess 候选 #1: 200.000 MHz ----------
[OUT #1] N=2048, Nr=2048, Fs_hat=409.600 MHz (基于 Fs_guess=200.000 MHz)

---------- Fs_guess 候选 #2: 240.000 MHz ----------
[OUT #2] N=4096, Nr=7004, Fs_hat=239.574 MHz (基于 Fs_guess=240.000 MHz)

---------- Fs_guess 候选 #3: 250.000 MHz ----------
[OUT #3] N=2048, Nr=3354, Fs_hat=250.000 MHz (基于 Fs_guess=250.000 MHz)

---------- Fs_guess 候选 #4: 300.000 MHz ----------
[OUT #4] N=1024, Nr=1398, Fs_hat=300.000 MHz (基于 Fs_guess=300.000 MHz)

========== 所有 Fs_guess 候选结果汇总 ==========
[候选#1] Fs_guess=200.000 MHz → N=2048, Fs_hat=409.600 MHz, dyn=12.3 dB
[候选#2] Fs_guess=240.000 MHz → N=4096, Fs_hat=239.574 MHz, dyn=15.8 dB
[候选#3] Fs_guess=250.000 MHz → N=2048, Fs_hat=250.000 MHz, dyn=11.2 dB
[候选#4] Fs_guess=300.000 MHz → N=1024, Fs_hat=300.000 MHz, dyn=9.5 dB

[最佳] 选择候选 #2 (Fs_guess=240.000 MHz, dyn=15.8 dB 最高)
```

## 查看所有候选详情

```matlab
% 遍历所有候选结果
for i = 1:numel(R.all_Fs_candidates)
    res = R.all_Fs_candidates{i};
    fprintf('候选 #%d:\n', i);
    fprintf('  Fs_guess  = %.3f MHz\n', res.Fs_guess/1e6);
    fprintf('  N_hat     = %d\n', res.N_hat);
    fprintf('  Fs_hat    = %.6f MHz\n', res.Fs_hat/1e6);
    fprintf('  动态范围  = %.2f dB\n', res.diagN.dyn_dB_best);
end
```

## 建议的候选值设置

根据不同应用场景，可以设置不同的候选值范围：

### Starlink 卫星信号
```matlab
opts.Fs_guess = [200e6, 240e6, 250e6];  % 200-250 MHz 范围
```

### 宽范围搜索
```matlab
opts.Fs_guess = [100e6, 150e6, 200e6, 240e6, 250e6, 300e6];
```

### 精细搜索（围绕已知值）
```matlab
opts.Fs_guess = (230:5:250)*1e6;  % 230-250 MHz，步进 5 MHz
```

## 性能注意事项

- 每个候选值都会完整运行一次估计算法
- 候选值越多，总运行时间越长
- 建议从少量候选开始（3-5个），根据结果再调整
- 使用 `opts.Verbose = false` 可减少输出，加快运行

## 参数调优

如果所有候选结果都不理想，可以调整以下参数：

```matlab
opts.p = 0.12;       % Sr 相对宽度，增大可扩展搜索范围
opts.nu_dB = 10;     % 验证门限，降低可放宽要求
opts.M_N = [];       % R0 计算样本数，可手动指定
```
