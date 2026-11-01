% ============================================================
% run_fig7.m — 论文 Figure 7（一维分岔图）的生成入口
%
% 这张图在论文里的位置：§4.4，Fig. 7（fig_bifurcation.pdf）
% 作用：沿 Gs 与 Ge 两条轴细扫，把"期"读成参数的分岔结构——主频在哪里跳、
%       幅度在哪里跳，与 v2 判据换色的位置是否重合。判据是读法，分岔是结构。
%
% 输入：无
% 输出：<项目根>/figures/fig_bifurcation.pdf
%
% 耗时：首次约 4 min（225 个参数点 × 70 s）；已有 results/scan_bifurcation.mat 时只画图。
% 依赖：tests/scan_bifurcation.m → plotting/fig7_bifurcation_paper.m
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

resfile = fullfile(root, 'tests', 'results', 'scan_bifurcation.mat');
if ~exist(resfile, 'file')
    fprintf('未找到 %s，先跑 tests/scan_bifurcation.m\n', resfile);
    scan_bifurcation();
end

outpath = fullfile(figdir, 'fig_bifurcation.pdf');
fig7_bifurcation_paper(outpath);
fprintf('fig7 已导出：%s\n', outpath);
