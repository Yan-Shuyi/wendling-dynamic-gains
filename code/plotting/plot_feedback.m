function plot_feedback(SavePath)
% PLOT_FEEDBACK  反向边 PZ→SOZ 的强度扫描：回灌加速点火，但不改熄火。
%
%   plot_feedback()            % 只画
%   plot_feedback(SavePath)    % 导出 PDF
%
% 依赖
%   tests/results/feedback.mat        w = 0/0.1/0.2/0.4
%   tests/results/feedback_strong.mat w = 0/0.6/0.8/1.0
%
% 两个面板
%   (a) 起振时刻随 w 的变化：SOZ 与 PZ 都提前，w=1.0 时 SOZ 提前 17.5 s。
%   (b) 峰值 rms 随 w 的变化：单调上升（+29%），但事件数始终是 1——
%       回灌改的是"什么时候点着"，不是"怎么熄火"。
if nargin < 1, SavePath = []; end

here = fileparts(mfilename('fullpath'));
resdir = fullfile(fileparts(here), 'tests', 'results');
A = load(fullfile(resdir, 'feedback.mat'));
B = load(fullfile(resdir, 'feedback_strong.mat'));

set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultLineLineWidth', 0.8);

col = [0.75 0.20 0.20; 0.90 0.55 0.15];
names = {'SOZ','PZ'};

% 合并两批，按 w 排序去重
ws = []; T1 = []; RM = [];
for S = {A, B}
    M = S{1};
    f = fieldnames(M);
    for k = 1:numel(f)
        e = M.(f{k});
        w = e.w;
        t1 = nan(1,2); rm = nan(1,2);
        for i = 1:2
            if ~isempty(e.ev{i}), t1(i) = e.ev{i}(1,1); end
            xh = highpass(e.y(i,:)', 0.5, 256)';   % 与判据同口径：先去直流
            w1 = 256; nr = floor(numel(xh)/w1);
            rr = zeros(1,nr);
            for j = 1:nr
                rr(j) = sqrt(mean(xh((j-1)*w1+(1:w1)).^2));
            end
            rm(i) = max(rr);
        end
        ws(end+1) = w; T1(end+1,:) = t1; RM(end+1,:) = rm; %#ok<AGROW>
    end
end
[ws, ix] = unique(ws);
T1 = T1(ix,:); RM = RM(ix,:);

f = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 6.3 2.7]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on');
h = gobjects(1,2);
for i = 1:2
    h(i) = plot(ax, ws, T1(:,i), '-o', 'Color', col(i,:), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 3.5);
end
xlim(ax, [0 1.05]); ylim(ax, [150 195]);
xlabel(ax, 'W(PZ \rightarrow SOZ)'); ylabel(ax, 'onset (s)');
legend(ax, h, names, 'Location', 'southwest', 'Box', 'off', 'FontSize', 7);
title(ax, '(a) onset advances with feedback', 'FontWeight', 'normal', 'FontSize', 8);
grid(ax, 'on'); set(ax, 'GridAlpha', 0.15);

ax = nexttile(tl); hold(ax, 'on');
h = gobjects(1,2);
for i = 1:2
    h(i) = plot(ax, ws, RM(:,i), '-s', 'Color', col(i,:), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 3.5);
end
xlim(ax, [0 1.05]); ylim(ax, [2.3 5.0]);
xlabel(ax, 'W(PZ \rightarrow SOZ)'); ylabel(ax, 'peak RMS (mV)');
legend(ax, h, names, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
title(ax, '(b) amplitude, not duration', 'FontWeight', 'normal', 'FontSize', 8);
grid(ax, 'on'); set(ax, 'GridAlpha', 0.15);

if ~isempty(SavePath)
    set(f, 'PaperUnits', 'inches', 'PaperSize', [6.3 2.7], ...
        'PaperPosition', [0 0 6.3 2.7], 'PaperPositionMode', 'manual');
    print(f, SavePath, '-dpdf', '-r600');
    fprintf('已导出 %s\n', SavePath);
end
end
