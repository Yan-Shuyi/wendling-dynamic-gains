% ============================================================
% run_fig5.m — 论文 Figure 5（事件对齐的分期统计）的生成入口
%
% 这张图在论文里的位置：§4.3，Fig. 5（fig_event_stats.pdf）
% 作用：回答"四个时期的增益统计量可不可分"。原做法在 10/30/50/70 s 四个固定时刻取样，
%       但 onset 在每一例里落在不同时刻，固定时刻取到的不是同一相位（D14 的结论）。
%       这里改成先按 onset 对齐、再在相对时刻取样，并把两种取法并列画出来。
%
% 输入：无
% 输出：<项目根>/figures/fig_event_stats.pdf
%
% 耗时：数据约 15 min（60 例 × 3 个触发幅度 × 200 s 仿真 + 20 例逐窗分类）；
%       已有 results/event_stats.mat 时只画图，几秒。
% 依赖：tests/event_stats.m → plotting/fig5_event_paper.m
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

resfile = fullfile(root, 'tests', 'results', 'event_stats.mat');
if ~exist(resfile, 'file')
    fprintf('未找到 %s，先跑 tests/event_stats.m\n', resfile);
    event_stats(60, 20);
end

outpath = fullfile(figdir, 'fig_event_stats.pdf');
fig5_event_paper(outpath);
fprintf('fig5 已导出：%s\n', outpath);
