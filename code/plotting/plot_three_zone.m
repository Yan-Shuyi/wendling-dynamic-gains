function plot_three_zone(SavePath)
% PLOT_THREE_ZONE  三区网络第一次自发启动：跨区次序图。
%
%   plot_three_zone()            % 只画
%   plot_three_zone(SavePath)    % 导出 PDF
%
% 依赖
%   tests/results/three_zone.mat（由 tests/run_three_zone.m 生成）
%
% 四个面板
%   (a) 三区 1 s 窗 rms 全时程：SOZ 先起，PZ 秒级跟上，NIZ 不动。
%   (b) 三区 G_e：启动器把 SOZ 的 G_e 从 3.25 缓慢推到 4.25；PZ 的 G_e 在被打响
%       之后同样爬到 8.5 —— 这是"募集"与"被动跟随"的分界。
%   (c) 原始 EEG：放电段三节点叠放（按各自 p99 幅度自适应偏移）。
%   (d) 170–340 s 的 G_s 放大：崩塌本身就是事件，三区的崩塌时刻差几秒。
if nargin < 1, SavePath = []; end

here = fileparts(mfilename('fullpath'));
resdir = fullfile(fileparts(here), 'tests', 'results');
L = load(fullfile(resdir, 'three_zone.mat'));

set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultLineLineWidth', 0.8);

col = [0.75 0.20 0.20;    % SOZ
       0.90 0.55 0.15;    % PZ
       0.20 0.42 0.70];   % NIZ
names = {'SOZ','PZ','NIZ'};

f = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 6.3 4.8]);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

D = L.net;  T = D.T;  fs = 256;
te = arrayfun(@(i) t_elec(D.yh(i,:), fs), 1:3);

%% (a) rms
ax = nexttile(tl); hold(ax, 'on');
h = gobjects(1,3);
for i = 1:3
    h(i) = plot(ax, D.tr, D.rr(i,:), 'Color', col(i,:), 'LineWidth', 0.9);
    if ~isnan(te(i))
        xline(ax, te(i), ':', 'Color', col(i,:), 'LineWidth', 0.7);
    end
end
xlim(ax, [0 T]); ylim(ax, [0 6]);
xlabel(ax, 'time (s)'); ylabel(ax, 'RMS (mV, 1 s, >0.5 Hz)');
legend(ax, h, names, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, '(a) activity in three zones', 'FontWeight', 'normal', 'FontSize', 8);

text(ax, 0.03, 0.55, sprintf('onset: SOZ %.0f s, PZ %+.0f s, NIZ none', ...
    te(1), te(2)-te(1)), 'Units', 'normalized', 'FontSize', 7);

%% (b) G_e
ax = nexttile(tl); hold(ax, 'on');
h = gobjects(1,3);
for i = 1:3
    h(i) = plot(ax, D.t, D.dy(3*(i-1)+1, :), 'Color', col(i,:), 'LineWidth', 0.9);
end
yline(ax, 4.0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.6);
xlim(ax, [0 T]); ylim(ax, [2.9 9]);
xlabel(ax, 'time (s)'); ylabel(ax, 'G_e (mV)');
legend(ax, h, names, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, '(b) excitatory gain', 'FontWeight', 'normal', 'FontSize', 8);

%% (c) raw EEG
ax = nexttile(tl); hold(ax, 'on');
[~, ir] = max(D.rr(1,:)); t0 = max(0, D.tr(ir) - 1.5);
w = round(3*fs); k1 = round(t0*fs) + 1; k2 = min(size(D.y,2), k1 + w - 1);
tt = (0:numel(k1:k2)-1)/fs;
% 先减窗内均值：放电期有一个很大的去极化偏置（21 mV 水平），
% 不去掉就只剩一条抬高的平台，看不见 3–4 Hz 的振荡。
seg = zeros(3, numel(k1:k2));
for i = 1:3
    v = D.yh(i, k1:k2);
    seg(i,:) = v - mean(v);
end
amp = zeros(1,3);
for i = 1:3
    amp(i) = prctile(abs(seg(i,:)), 99);
end
lv = zeros(1,3);
lv(3) = 0;
lv(2) = 1.8*max(amp(3), 0.45);
lv(1) = lv(2) + 1.8*max(amp(2), 2);
for i = 1:3
    plot(ax, tt, seg(i,:) + lv(i), 'Color', col(i,:), 'LineWidth', 0.5);
end
xlim(ax, [0 tt(end)]);
ylim(ax, [-1.5*max(amp(3),0.45), lv(1) + 1.5*amp(1)]);
yticklabels(ax, fliplr(names));   % 刻度自下而上，SOZ 在最上面（偏移最大）
xlabel(ax, sprintf('time in %.1f–%.1f s window (s)', t0, t0+3));
title(ax, '(c) raw EEG, 3 s of the discharge', 'FontWeight', 'normal', 'FontSize', 8);

%% (d) G_s 放大
ax = nexttile(tl); hold(ax, 'on');
h = gobjects(1,3);
for i = 1:3
    h(i) = plot(ax, D.t, D.dy(3*(i-1)+2, :), 'Color', col(i,:), 'LineWidth', 0.9);
end
xlim(ax, [170 340]); ylim(ax, [0 36]);
xlabel(ax, 'time (s)'); ylabel(ax, 'G_s (mV)');
legend(ax, h, names, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
title(ax, '(d) slow inhibitory gain', 'FontWeight', 'normal', 'FontSize', 8);

if ~isempty(SavePath)
    set(f, 'PaperUnits', 'inches', 'PaperSize', [6.3 4.8], ...
        'PaperPosition', [0 0 6.3 4.8], 'PaperPositionMode', 'manual');
    print(f, SavePath, '-dpdf', '-r600');
    fprintf('已导出 %s\n', SavePath);
end
end

function t0 = t_elec(x, fs)
% 0.25 s 窗 rms，首次持续 2 s 越过 max(4×基线, 0.5 mV)
w = round(0.25*fs); nr = floor(numel(x)/w);
rr = zeros(1, nr);
for j = 1:nr
    seg = x((j-1)*w+(1:w));
    rr(j) = sqrt(mean(seg.^2));
end
tr = (0:nr-1)*0.25 + 0.125;
base = median(rr(tr <= 30));
above = rr > max(4*base, 0.5);
t0 = NaN;
for j = 1:nr-8
    if all(above(j:j+7)), t0 = tr(j)-0.125; return; end
end
end
