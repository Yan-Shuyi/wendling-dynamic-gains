% ============================================================
% run_fig1.m — 论文 Figure 3（增益扫描图）的生成入口
%
% 这张图在论文里的位置：§4.1，Fig. 3（fig_gain_sweep.pdf）
% 作用：7 个互不耦合的单节点仿真，每个节点取一个慢抑制增益 G_s（5:5:35），
%       其余参数固定，用来展示"仅改变抑制强度就能把系统推向不同动力学状态"。
%       上排是时域波形，下排是对应 PSD，背景色块标出 δ/θ/α/β/γ 频段。
%
% 输入：无（参数在下面"参数"区直接改）
% 输出：<项目根>/figures/fig_gain_sweep.pdf
%
% 耗时：约 21 s（7 节点 × 200 s 仿真）
% 依赖：default_config → simulate_network → fig1_paper
%
% 注：这是 7 次独立仿真（sc_matrix 全零），不是 7 节点耦合网络。
%     论文正文必须写成 seven independent single-node simulations，
%     写成 "7-node coupled network" 是早期错误。
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

%% 参数
cfg = default_config();
cfg.node_count    = 7;
cfg.A             = 4.25;
cfg.B             = 5:5:35;     % B 即 G_s：逐节点取不同慢抑制增益
cfg.G             = 10;         % G 即 G_f
cfg.model_dynamic = 0;          % 静态增益（不启用 ESF 迟滞律）
cfg.esf_enable    = false;      % 不做病理触发
cfg.resample_dy   = false;      % 不需要 E/S/F 轨迹

%% 仿真
fprintf('fig1: %d 个节点，G_s = %s\n', cfg.node_count, mat2str(cfg.B));
tic;
outputs_fs = simulate_network(cfg);
fprintf('fig1 仿真完成，耗时 %.1f s\n', toc);

%% 出图
outpath = fullfile(figdir, 'fig_gain_sweep.pdf');
fig1_paper(cfg.node_count, cfg.simulation_time, cfg.eeg_fs, ...
           outputs_fs, cfg.A, cfg.B, cfg.G, ...
           't_skip', 1, 't_view', 5, 'f_max', 40, ...
           'SavePath', outpath);

fprintf('fig1 已导出：%s\n', outpath);
