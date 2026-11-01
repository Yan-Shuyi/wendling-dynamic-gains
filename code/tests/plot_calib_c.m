function plot_calib_c(R)
% PLOT_CALIB_C  画 C 标定曲线。R 来自 tests/results/calib_c_sweep.mat
%
%   plot_calib_c()     从 mat 读
%   plot_calib_c(R)    直接给结构体
%
% 四个面板：rms、rms 比、主频、PLV，横轴 C（含 C=0 点），竖虚线标 C=20。

if nargin < 1 || isempty(R)
    root = fileparts(fileparts(mfilename('fullpath')));
    S = load(fullfile(root, 'tests', 'results', 'calib_c_sweep.mat'));
    R = S.R;
else
    root = fileparts(fileparts(mfilename('fullpath')));
end
Cv = R.Cv;
names = R.recv(:, 1);
cols = [0.00 0.45 0.70; 0.85 0.33 0.10];
Cref = 15;

f = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 6.3 4.0], 'Visible', 'on');
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

items = {'rms', 'ratio', 'fr', 'plv'};
tits = {'(a) receiver rms (mV)', '(b) rms / uncoupled', ...
        '(c) dominant frequency (Hz)', '(d) phase-locking with driver'};
ylabs = {'', 'ratio', 'Hz', 'PLV'};
for k = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    for m = 1:2
        nm = names{m};
        switch items{k}
            case 'rms',   yy = R.(nm).rms;
            case 'ratio', yy = R.(nm).rms ./ R.(nm).rms(1);
            case 'fr',    yy = R.(nm).fr;
            case 'plv',   yy = R.(nm).plv;
        end
        yy = yy(:)';
        plot(ax, Cv, yy, '-o', 'Color', cols(m,:), 'MarkerFaceColor', 'w', ...
            'LineWidth', 1.3, 'MarkerSize', 4.5, 'DisplayName', nm);
    end
    xline(ax, Cref, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xlabel(ax, 'coupling constant C');
    title(ax, tits{k}, 'FontWeight', 'normal', 'FontSize', 9);
    grid(ax, 'on'); box(ax, 'off');
    ax.Toolbar.Visible = 'off';
    xlim(ax, [0 120]); xticks(ax, [0 20 40 60 80 100 120]);
    if ~isempty(ylabs{k}), ylabel(ax, ylabs{k}); end
    if k == 2, set(ax, 'YScale', 'log'); end
    if k == 1, legend(ax, 'Location', 'northwest', 'Box', 'off'); end
    if k == 3, ylim(ax, [0 12]); end
    if k == 4, ylim(ax, [0.4 1.0]); end
end

figdir = fullfile(root, '..', 'figures');
exportgraphics(f, fullfile(figdir, 'calib_coupling_c.png'), 'Resolution', 200);
close(f);
fprintf('已写 figures/calib_coupling_c.png\n');
end
