% ============================================================
% run_fig2.m — 论文 Figure 4（单例演化图）的生成入口
%
% 这张图在论文里的位置：§4.2，Fig. 4（fig_case_evolution.pdf）
% 作用：单节点、单次仿真，t = 10 s 时把 G_e 从 3.25 跳到 4.25 mV，
%       即人为触发一次 ESF。图分三层：
%         上：EEG 波形（静息 → 高幅慢节律放电 → 崩塌后抑制 → 恢复）
%         中：v2 判据（§3.5、Fig. 2）给出的四期色条
%         下：G_e / G_s / G_f 三条轨迹
%       竖直虚线两条：t = 10 s（扰动）与 G_s 崩塌时刻（onset）。
%
% 输入：无（参数在下面"参数"区直接改）
% 输出：<项目根>/figures/fig_case_evolution.pdf
%
% 耗时：约 5 s
% 依赖：default_config → simulate_network → classify_params_v2 /
%       make_reference / classify_timeseries_v2 → plot_case_v2
%
% 注：2026-10-31 起改用 v2 判据。旧版走 classify_eeg_timeseries（5 s 窗、六类
%     振幅分区、Type 1–6），与论文 §3.5 的频带判据不是一套，已停用。
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

%% 参数
cfg = default_config();
cfg.noise_segments = [0.3 0.5 0.2];   % 历史设置：本图结果基于三段噪声
cfg.node_count  = 1;
cfg.A           = 3.25;             % 初始 G_e
cfg.B           = 22;               % G_s
cfg.G           = 10;               % G_f
cfg.esf_enable  = true;             % 启用动态增益（ESF 迟滞律）
cfg.esf_time    = 10;               % 第 10 s 触发
cfg.esf_params  = [4.25, 22, 10];   % 触发后跳到的 (G_e, G_s, G_f)
cfg.esf_jitter  = 0;                % 单例固定触发，不加随机抖动
cfg.resample_dy = true;

%% 仿真
fprintf('fig2: 单节点，%.0f s，t = %.0f s 触发 ESF\n', cfg.simulation_time, cfg.esf_time);
tic;
[out_fs, dy_fs, cfg] = simulate_network(cfg);
fprintf('fig2 仿真完成，耗时 %.1f s\n', toc);

x  = out_fs(1, :);
Fs = cfg.eeg_fs;
t  = (0:numel(x)-1) / Fs;

%% 分期：参考段取扰动前 10 s（此时增益仍在静息值）
p   = classify_params_v2();
ref = make_reference(x(1:round(cfg.esf_time*Fs)), Fs, p);
[~, labs, ~] = classify_timeseries_v2(x, Fs, p, ref);

%% onset：G_s 首次跌破中点（与 tests/event_stats.m 同一约定）
thr  = 12.75;
tg   = (0:size(dy_fs,2)-1) / cfg.dy_fs;
k    = find(dy_fs(2, :) < thr, 1, 'first');
t_on = tg(k);

fprintf('onset（G_s < %.2f）= %.1f s\n', thr, t_on);
fprintf('出现的标签：%s\n', strjoin(unique(labs), ', '));

%% 出图
tc = (0:2:(cfg.simulation_time-10)) + 5;      % 每窗中心（10 s 窗 / 2 s 步）
outpath = fullfile(figdir, 'fig_case_evolution.pdf');
plot_case_v2(x, t, labs, tc, dy_fs(1, :), dy_fs(2, :), dy_fs(3, :), Fs, p, ...
             [cfg.esf_time, t_on], outpath);
fprintf('fig2 已导出：%s\n', outpath);
