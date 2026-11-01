% ============================================================
% run_all.m — 一键重跑论文全部仿真图
%
% 依次执行，图统一输出到 <项目根>/figures/：
%   run_fig2          fig_case_evolution.pdf   （论文 Fig. 4，单例周期）
%   run_fig5          fig_event_stats.pdf      （论文 Fig. 5，事件对齐统计）
%   run_fig7          fig_bifurcation.pdf      （论文 Fig. 6，一维分岔图）
%   run_fig3_phase    fig_phase_BG.pdf         （论文 Fig. 3，(Gs,Gf) 参数平面）
%
% 论文 Fig. 1（模型结构）与 Fig. 2（判定链）不是 MATLAB 生成的：它们是 TikZ 源码
% （figures/fig_model_structure.tex、fig_classifier_tree.tex），改完在该目录下用
% pdflatex 编译即可。
%
% 数据依赖（首次跑会各花几分钟，之后直接读 .mat）：
%   tests/event_stats.m   → results/event_stats.mat  （60 例 × 3 个触发幅度）
%   tests/event_labels.m  → results/event_labels.mat （逐窗标签）
%   tests/scan_bifurcation.m('S1'/'S2') → results/bifur_S1.mat、bifur_S2.mat
%   tests/scan_bg.m       → results/scan_bg_A*.mat   （三个 A 的参数平面）
%
% 总耗时（无缓存）约 40 分钟，主要在仿真；建议在 MATLAB 里直接跑
%
%   >> run_all
%
% 而不要求从某个特定目录启动——脚本会自己 cd 到 code/。
%
% 注：run_fig1 / run_fig3（旧的 Gs 扫描图与旧箱线图）已不再进入论文，
% 保留脚本仅用于追溯；论文 Fig. 5 已由事件对齐版替代。

here = fileparts(mfilename('fullpath'));
cd(here);

t_all = tic;

fprintf('\n===== [1/4] fig2（论文 Fig. 4，单例周期）=====\n');
run('run_fig2.m');

fprintf('\n===== [2/4] fig5（论文 Fig. 5，事件对齐统计）=====\n');
run('run_fig5.m');

fprintf('\n===== [3/4] fig7（论文 Fig. 6，分岔图）=====\n');
run('run_fig7.m');

fprintf('\n===== [4/4] fig3_phase（论文 Fig. 3，参数平面）=====\n');
run('run_fig3_phase.m');

fprintf('\n全部完成，总耗时 %.1f s。输出目录：%s\n', ...
        toc(t_all), fullfile(here, '..', 'figures'));
