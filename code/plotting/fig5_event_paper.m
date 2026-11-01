function fig5_event_paper(SavePath)
% FIG5_EVENT_PAPER  论文 Fig.5：分期统计必须锚在 onset 上，不能锚在钟上。
%
%   fig5_event_paper()            % 只画
%   fig5_event_paper(SavePath)    % 导出 PDF
%
% 输入
%   SavePath   导出路径；省略则只显示
%
% 依赖
%   tests/results/event_stats.mat（由 tests/event_stats.m 生成）
%
% 四个面板
%   (a) 60 例的 Gs 轨迹对齐到 onset：灰线是各例，粗线是中位数，阴影是 IQR 带。
%       对齐之后崩塌落在同一竖线上——这是"onset 是事件不是时刻"的直接图像。
%   (b) 对齐取样：四个相位 × 三个触发幅度的 Gs 箱线图。三色几乎重合，说明
%       相位统计量不随触发幅度变。
%   (c) 固定时刻取样：同样四个"相位"（9/30/50/70 s）。同一个箱子里三色分开，
%       因为不同幅度的 onset 落在不同时刻，固定时刻取到的不是同一相位。
%   (d) 相对 onset 的标签构成（v2 判据，10 s 窗 / 2 s 步）。发作起始处由
%       interictal 切换到 ONSET，增益崩塌后回到慢节律/背景。
if nargin < 1, SavePath = []; end

here = fileparts(mfilename('fullpath'));
resdir = fullfile(fileparts(here), 'tests', 'results');
D  = load(fullfile(resdir, 'event_stats.mat'));
DL = load(fullfile(resdir, 'event_labels.mat'));   % (d) 用逐窗标签（参考段更干净）

set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultLineLineWidth', 0.8);

cmp = [0.86 0.86 0.86;   % intICTAL
       0.96 0.78 0.45;   % preONSET
       0.85 0.45 0.75;   % ONSET
       0.25 0.50 0.80;   % ICTAL
       0.55 0.35 0.65;   % postICTAL（默认判据不下发，保留配色以便扩展）
       0.60 0.60 0.60];  % uncertain
acol = [0.20 0.35 0.65; 0.85 0.45 0.15; 0.30 0.60 0.35];   % 三个触发幅度

f = figure('Visible','off', 'Color','w', 'Units','inches', ...
    'Position',[1 1 6.3 5.0], 'PaperUnits','inches', ...
    'PaperSize',[6.3 5.0], 'PaperPosition',[0 0 6.3 5.0], ...
    'PaperPositionMode','manual');
tl = tiledlayout(2, 2, 'Padding','compact', 'TileSpacing','compact');

%% ---------- (a) 对齐后的 Gs 轨迹 ----------
a0 = find(abs(D.ges - 4.25) < 1e-9, 1);
if isempty(a0), a0 = 1; end
S  = D.traj{a0};
tg = D.tg;
ax = nexttile(tl); hold(ax, 'on');
med = median(S, 1, 'omitnan');
lo  = prctile(S, 25, 1);  hi = prctile(S, 75, 1);
fill(ax, [tg fliplr(tg)], [lo fliplr(hi)], [0.85 0.88 0.93], ...
    'EdgeColor','none');
for i = 1:size(S, 1)
    plot(ax, tg, S(i, :), 'Color', [0.72 0.74 0.78], 'LineWidth', 0.4);
end
plot(ax, tg, med, 'k', 'LineWidth', 1.4);
xl = [tg(1) tg(end)];
plot(ax, [0 0], [0 32], 'k--', 'LineWidth', 0.7);
text(ax, 1.5, 30, 'collapse', 'FontSize', 7);
ylim(ax, [0 32]); xlim(ax, xl);
xlabel(ax, 'time relative to the collapse (s)');
ylabel(ax, 'G_s (mV)');
title(ax, sprintf('(a) G_s aligned to the collapse, %d runs, G_e = %.2f', size(S, 1), D.ges(a0)), ...
    'FontWeight','normal');

%% ---------- (b) 对齐取样：相位 × 幅度 ----------
ax = nexttile(tl); hold(ax, 'on');
draw_box(ax, D.al_all, D.ges, acol);
set(ax, 'XTick', 1:4, 'XTickLabel', {'interictal','pre-collapse','collapse','post-collapse'}, ...
    'XTickLabelRotation', 15, 'FontSize', 7);
ylabel(ax, 'G_s (mV)');
title(ax, '(b) sampled at fixed phase relative to the collapse', 'FontWeight','normal');
ylim(ax, [0 32]);

%% ---------- (c) 固定时刻取样 ----------
ax = nexttile(tl); hold(ax, 'on');
draw_box(ax, D.fx_all, D.ges, acol);
set(ax, 'XTick', 1:4, 'XTickLabel', {'9 s','30 s','50 s','70 s'});
ylabel(ax, 'G_s (mV)');
xlabel(ax, 'absolute time (s)');
title(ax, '(c) sampled at fixed times', 'FontWeight','normal');
ylim(ax, [0 32]);

lg = legend(ax, arrayfun(@(g) sprintf('G_e = %.2f', g), D.ges, 'UniformOutput', false), ...
    'Orientation','horizontal', 'Box','off', 'FontSize', 7);
lg.Layout.Tile = 'south';

%% ---------- (d) 标签构成 ----------
ax = nexttile(tl); hold(ax, 'on');
% bin 宽度不等（20/10/10/10/20/40/80 s），所以按真实宽度画，横轴才是时间
cnt = DL.cnt;  cnt(isnan(cnt)) = 0;
frac = cnt ./ max(sum(cnt, 2), 1);
cum = zeros(size(frac,1), 1);
for q = 1:size(frac, 2)
    lo = cum;  hi = cum + frac(:, q);
    cum = hi;
    for j = 1:size(frac, 1)
        if frac(j, q) <= 0, continue; end
        patch(ax, [DL.edges(j) DL.edges(j+1) DL.edges(j+1) DL.edges(j)], ...
            [lo(j) lo(j) hi(j) hi(j)], cmp(q, :), 'EdgeColor', 'none');
    end
end
set(ax, 'XLim', [DL.edges(1) DL.edges(end)], 'YLim', [0 1]);
xlabel(ax, 'window centre relative to the collapse (s)');
ylabel(ax, 'fraction of windows');
title(ax, '(d) label composition (10-s windows)', 'FontWeight','normal');
% 图例要用显式句柄：patch 的叠放次序与标签次序不一致，默认图例会张冠李戴
hs = gobjects(1, numel(DL.labs));
for q = 1:numel(DL.labs)
    hs(q) = patch(ax, NaN, NaN, cmp(q, :), 'EdgeColor', 'none');
end
lg2 = legend(ax, hs, DL.labs, 'Orientation','horizontal', 'Box','off', 'FontSize', 6.5);
lg2.Layout.Tile = 'south';
lg2.NumColumns = 3;

if ~isempty(SavePath)
    print(f, SavePath, '-dpdf');
    fprintf('图已写：%s\n', SavePath);
end
end

% ----------------------------------------------------------------------
function draw_box(ax, A, ges, acol)
% A{a}{gain}：第 a 个幅度、第 g 个增益的 相位 × 例 矩阵。这里只画 Gs（g=2）。
w = 0.24;
for a = 1:numel(ges)
    M = A{a, 2};                      % Gs
    for q = 1:size(M, 1)
        v = M(q, ~isnan(M(q, :)));
        if isempty(v), continue; end
        xc = q + (a - (numel(ges)+1)/2) * w;
        boxchart(ax, xc*ones(size(v)), v, 'BoxWidth', w*0.85, ...
            'BoxFaceColor', acol(a, :), 'MarkerStyle','none', ...
            'WhiskerLineColor', acol(a, :), 'LineWidth', 0.6);
    end
end
xlim(ax, [0.5 size(A{1,2},1)+0.5]);
end
