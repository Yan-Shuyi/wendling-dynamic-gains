function fig7_bifurcation_paper(SavePath)
% FIG7_BIFURCATION_PAPER  论文 Fig.7：把"期"读成分岔结构。
%
%   fig7_bifurcation_paper()            % 只画
%   fig7_bifurcation_paper(SavePath)    % 导出 PDF
%
% 依赖
%   tests/results/scan_bifurcation.mat（由 tests/scan_bifurcation.m 生成）
%
% 四个面板
%   (a) 主频随慢抑制增益 Gs 的变化（Gf = 0 / 10 / 25 三条）。每个点是某个参数点上
%       一个 10 s 窗的谱峰；同一参数的竖直散布就是窗间抖动。
%   (b) 同一扫描的 rms。快活动支电压低、慢节律支电压高，两支在幅度上也分得开。
%   (c)(d) 主频与 rms 随兴奋性增益 Ge 的变化（两组抑制参数）：慢节律只在 Ge 超过
%       某个值之后出现，这正是 §4.2 那句"Ge 决定有没有慢节律"的图形版本。
%
% 读法
%   两个面板的横轴都是参数，不是时间；颜色是 v2 判据给的标签。判据的阈值线在图上
%   完全不可见，但颜色的切换位置与主频跳变的位置重合——这就是"判据只是读法"。
if nargin < 1, SavePath = []; end

here = fileparts(mfilename('fullpath'));
resdir = fullfile(fileparts(here), 'tests', 'results');
D1 = load(fullfile(resdir, 'bifur_S1.mat'));
D2 = load(fullfile(resdir, 'bifur_S2.mat'));

set(groot, ...
    'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize', 8, 'DefaultLineLineWidth', 0.9);

cmap = [0.86 0.86 0.86;   % intICTAL
        0.96 0.78 0.45;   % preONSET
        0.85 0.45 0.75;   % ONSET
        0.25 0.50 0.80;   % ICTAL
        0.60 0.60 0.60];  % uncertain
clab = {'intICTAL','preONSET','ONSET','ICTAL','uncertain'};
lcol = [0.20 0.35 0.65; 0.85 0.45 0.15; 0.30 0.60 0.35];

f = figure('Visible','off', 'Color','w', 'Units','inches', ...
    'Position',[1 1 6.3 3.7], 'PaperUnits','inches', ...
    'PaperSize',[6.3 3.7], 'PaperPosition',[0 0 6.3 3.7], ...
    'PaperPositionMode','manual');
tl = tiledlayout(2, 2, 'Padding','compact', 'TileSpacing','compact');

%% ---------- 上排：Gs 扫描 ----------
ax = nexttile(tl); hold(ax,'on');
scatter_by_label(ax, D1.F, D1.S, cmap, clab, 'f', [0 32]);
ylabel(ax, 'dominant frequency (Hz)');
xlabel(ax, 'G_s (mV)');
title(ax, {'(a) G_s sweep, G_e = 4.60 mV'; 'solid: G_f = 10,  dashed: G_f = 30'}, 'FontWeight','normal');
yline(ax, 15, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);

ax = nexttile(tl); hold(ax,'on');
scatter_by_label(ax, D1.F, D1.S, cmap, clab, 'rms', [0 7.5]);
ylabel(ax, 'rms (mV)');
xlabel(ax, 'G_s (mV)');
title(ax, '(b) same sweep, amplitude', 'FontWeight','normal');

%% ---------- 下排：Ge 扫描 ----------
ax = nexttile(tl); hold(ax,'on');
scatter_by_label(ax, D2.F, D2.S, cmap, clab, 'f', [0 13]);
ylabel(ax, 'dominant frequency (Hz)');
xlabel(ax, 'G_e (mV)');
title(ax, {'(c) G_e sweep'; 'solid: (G_s,G_f)=(25,10),  dashed: (10,10)'}, 'FontWeight','normal');
yline(ax, 15, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);

ax = nexttile(tl); hold(ax,'on');
scatter_by_label(ax, D2.F, D2.S, cmap, clab, 'rms', [0 7.5]);
ylabel(ax, 'rms (mV)');
xlabel(ax, 'G_e (mV)');
title(ax, '(d) same sweep, amplitude', 'FontWeight','normal');

% 图例：类别（颜色）与曲线（Gf / 抑制参数）分开说明
h = gobjects(numel(clab), 1);
for q = 1:numel(clab)
    h(q) = scatter(ax, NaN, NaN, 18, cmap(q,:), 'filled');
end
lg = legend(ax, h, clab, 'Orientation','horizontal', 'Box','off', 'FontSize', 7);
lg.Layout.Tile = 'south';

if ~isempty(SavePath)
    print(f, SavePath, '-dpdf');
    fprintf('图已写：%s\n', SavePath);
end
end

% ----------------------------------------------------------------------
function scatter_by_label(ax, F, S, cmap, clab, field, ylims)
% F 是 extract_all 的结果：f/rms/rf/rs/lab 是 nx × ns × nw
[nx, ns, nw] = size(F.f);
for j = 1:ns
    xs = []; ys = []; cs = [];
    for i = 1:nx
        for w = 1:nw
            v = F.(field)(i, j, w);
            if isnan(v), continue; end
            xs(end+1) = S.x(i);            %#ok<SAGROW>
            ys(end+1) = v;                 %#ok<SAGROW>
            q = find(strcmp(clab, F.lab{i,j,w}), 1);
            if isempty(q), q = numel(clab); end
            cs(end+1, :) = cmap(q, :);     %#ok<SAGROW>
        end
    end
    if isempty(xs), continue; end
    scatter(ax, xs, ys, 5, cs, 'filled', 'MarkerFaceAlpha', 0.55);
    % 中位数折线由调用方统一画（需要按曲线区分线型），这里只放散点
end
xlim(ax, [min(S.x) max(S.x)]);
ylim(ax, ylims);
% 曲线用线型区分（标题里注明），不再逐条标字：两条背景支会重合，标字必然叠在一起
ls = {'-','--'};
for j = 1:min(numel(S.lbl), 2)
    yv = squeeze(median(F.(field)(:, j, :), 3, 'omitnan'));
    plot(ax, S.x(:), yv(:), ls{j}, 'Color', [0.10 0.10 0.10], 'LineWidth', 1.1);
end
end
