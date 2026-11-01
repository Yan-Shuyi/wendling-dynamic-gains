% ============================================================
% run_fig3.m — 论文 Figure 5（分期箱线图）的生成入口
%
% 这张图在论文里的位置：§4.3，Fig. 5（fig_esf_boxplot.pdf）
% 作用：100 次互相解耦（sc_matrix 全零）的单节点仿真，每次触发增益带 ±0.05 的
%       随机抖动，相当于 100 次独立实验。把每次实验在 10 / 30 / 50 / 70 s 四个
%       时间点（对应 interICTAL / preONSET / ONSET / ICTAL）的
%       G_e（EXC）、G_s（SDI）、G_f（FSI）抽出来画箱线图，
%       用统计分布而不是单条轨迹说明"四个时期增益是可分的"。
%
% 输入：无（参数在下面"参数"区直接改）
% 输出：<项目根>/figures/fig_esf_boxplot.pdf
%
% 耗时：约 5 分 10 秒（100 节点 × 200 s 仿真 + 重采样，是三个脚本里最慢的）
% 依赖：default_config → simulate_network → fig3_paper
%
% 注：论文图注写的是 "100 independent simulation runs"，对应代码里
%     node_count = 100 且 sc_matrix 全零——即 100 次相互独立的实验，
%     不是 100 节点耦合网络。
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

%% 参数
cfg = default_config();
cfg.noise_segments = [0.3 0.5 0.2];   % 历史设置：本图结果基于三段噪声
cfg.node_count = 100;           % 100 次独立实验
cfg.A          = 3.25;          % 初始 G_e
cfg.B          = 22;            % G_s
cfg.G          = 10;            % G_f
cfg.esf_enable = true;
cfg.esf_time   = 10;            % 第 10 s 触发
cfg.esf_params = [4.25, 22, 10];% 触发后的目标增益
cfg.esf_jitter = 0.05;          % 每次实验的触发增益略有差异 → 制造样本差异

%% 仿真
fprintf('fig3: %d 次独立仿真，%.0f s\n', cfg.node_count, cfg.simulation_time);
tic;
[~, dy_fs, cfg] = simulate_network(cfg);
fprintf('fig3 仿真完成，耗时 %.1f s\n', toc);

%% 出图
outpath = fullfile(figdir, 'fig_esf_boxplot.pdf');
fig3_paper(dy_fs, cfg.dy_fs, [], [], outpath);

fprintf('fig3 已导出：%s\n', outpath);
