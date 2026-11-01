% ============================================================
% run_fig3_phase.m — 论文 Fig.3（(B,G) 参数平面分区图）的生成入口
%
% 这张图在论文里的位置：§4.1（替代原来"时间轴上分三相"的 Fig.3）
% 作用：把"期"画在参数轴上。三个子图是 A = 3.25 / 4.00 / 4.60 三个平面，
%       可以看出 (i) 快活动区固定在低 B、高 G 的角落，(ii) 慢节律区随 A 增大而扩张，
%       (iii) 两区之间隔着背景/不确定带。下排给出连续场与代表波形。
%
% 输入：无（数据直接读 code/tests/results/scan_bg_A*.mat）
% 输出：<项目根>/figures/fig_phase_BG.pdf
%
% 前置条件：先跑过参数扫描
%   addpath('C:/Users/26442/WorkBuddy/Claw/code'); setup_paths(); scan_bg([],[],3.25);
%   scan_bg([],[],4.00); scan_bg([],[],4.60)
%
% 依赖：fig3_phase_paper → results/scan_bg_<tag>.mat → simulate_network（画波形用）
% ============================================================

clear; clc; close all;
[root, figdir] = setup_paths();

outpath = fullfile(figdir, 'fig_phase_BG.pdf');
fig3_phase_paper({'A3p25','A4p00','A4p60'}, outpath);

fprintf('Fig.3（参数平面分区）已导出：%s\n', outpath);
