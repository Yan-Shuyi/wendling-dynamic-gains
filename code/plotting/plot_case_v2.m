function plot_case_v2(x, t, labs, t_center, Ge, Gs, Gf, Fs, p, marks, SavePath)
% PLOT_CASE_V2  论文 Fig.4：EEG / 四期色条 / 三条增益轨迹（三层共享时间轴）。
%
%   plot_case_v2(x, t, labs, t_center, Ge, Gs, Gf, Fs, p, marks, SavePath)
%
% 输入
%   x         整段 EEG（1 × N，mV）
%   t         对应时间轴（s）
%   labs      每窗标签（cell，来自 classify_timeseries_v2）
%   t_center  每窗中心时刻（s）
%   Ge/Gs/Gf  三条增益轨迹（1 × n，采样率 Fs）
%   Fs        EEG 采样率（Hz）
%   p         classify_params_v2()
%   marks     竖直参考线时刻（s），如 [esf_time, collapse_time]
%   SavePath  导出路径（省略则只画）
%
% 备注
%   色条按窗步长画成 2 s 宽的色块，所以色条的边界就是窗边界，读图不会误以为标签
%   有比 2 s 更细的时间分辨率。参考线用 HandleVisibility='off'，不进图例。
if nargin < 11, SavePath = []; end

set(groot, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize', 8);

cmap = [0.94 0.94 0.94;   % intICTAL
        0.96 0.78 0.45;   % preONSET
        0.85 0.45 0.75;   % ONSET
        0.25 0.50 0.80;   % ICTAL
        0.55 0.35 0.65;   % postICTAL
        0.60 0.60 0.60];  % uncertain

f = figure('Visible','off', 'Color','w', 'Units','inches', ...
    'Position',[1 1 6.3 4.0], 'PaperUnits','inches', ...
    'PaperSize',[6.3 4.0], 'PaperPosition',[0 0 6.3 4.0], ...
    'PaperPositionMode','manual');

ax1 = axes('Parent', f, 'Position', [0.11 0.585 0.85 0.375]);
ax2 = axes('Parent', f, 'Position', [0.11 0.495 0.85 0.070]);
ax3 = axes('Parent', f, 'Position', [0.11 0.115 0.85 0.290]);

%% 上：EEG（按列取 min/max 包络画成带子）
% 5 万个采样点压到 5 in 宽，直接画线在 PDF 里会糊成一块或干脆看不见；
% 按"每列 min/max"画成填充带，振幅包络与直流偏移都保留，且矢量文件很小。
hold(ax1, 'on');
ncol = 900;
idx = round(linspace(1, numel(x), ncol+1));
xe = zeros(1, ncol); lo = zeros(1, ncol); hi = zeros(1, ncol);
for c = 1:ncol
    k1 = idx(c); k2 = max(idx(c+1)-1, k1);
    xe(c) = mean(t(k1:k2));
    lo(c) = min(x(k1:k2));  hi(c) = max(x(k1:k2));
end
patch(ax1, [xe fliplr(xe)], [lo fliplr(hi)], [0.15 0.15 0.15], ...
    'EdgeColor', 'none');
ylabel(ax1, 'EEG (mV)');
yl = prctile(x, [0.2 99.8]);
ylim(ax1, yl + 0.08*[-1 1]*diff(yl));
set(ax1, 'XTickLabel', [], 'Box', 'on', 'XLim', [t(1) t(end)], 'TickDir', 'out');

%% 中：标签色条
hold(ax2, 'on');
nw = numel(labs);
for w = 1:nw
    li = find(strcmp(p.names, labs{w}), 1);
    if isempty(li), li = numel(p.names); end
    tc = t_center(w);
    patch(ax2, [tc-1 tc+1 tc+1 tc-1], [0 0 1 1], cmap(li, :), 'EdgeColor', 'none');
end
set(ax2, 'XLim', [t(1) t(end)], 'YLim', [0 1], 'YTick', [], ...
    'XTickLabel', [], 'Box', 'on');
ylabel(ax2, 'state', 'FontSize', 7);

%% 下：三条增益
hold(ax3, 'on');
tg = (0:numel(Ge)-1)/Fs;
h1 = plot(ax3, tg, Ge, '-', 'Color', [0.75 0.25 0.25], 'LineWidth', 0.9);
h2 = plot(ax3, tg, Gs, '-', 'Color', [0.20 0.40 0.75], 'LineWidth', 0.9);
h3 = plot(ax3, tg, Gf, '-', 'Color', [0.25 0.60 0.35], 'LineWidth', 0.9);
ylabel(ax3, 'synaptic gain (mV)');
xlabel(ax3, 'time (s)');
set(ax3, 'Box', 'on', 'XLim', [t(1) t(end)], 'TickDir', 'out');
% 曲线名直接标在右端（三条曲线在 200 s 处相距足够远），比图例省地方也更好读
xl3 = t(end);
for hh = {h1, 'G_e'; h2, 'G_s'; h3, 'G_f'}'
    yv = hh{1}.YData(end);
    text(ax3, xl3 - 0.01*diff(get(ax3,'XLim')), yv, hh{2}, ...
        'Color', hh{1}.Color, 'FontSize', 7, 'FontWeight','bold', ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle');
end

%% 参考线（扰动 / 崩塌），不参与图例
for m = marks(:)'
    plot(ax1, [m m], get(ax1,'YLim'), ':', 'Color', [0.30 0.30 0.30], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax2, [m m], [0 1], ':', 'Color', [0.30 0.30 0.30], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax3, [m m], get(ax3,'YLim'), ':', 'Color', [0.30 0.30 0.30], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
end

if ~isempty(SavePath)
    print(f, SavePath, '-dpdf');
    fprintf('图已写：%s\n', SavePath);
end
end
