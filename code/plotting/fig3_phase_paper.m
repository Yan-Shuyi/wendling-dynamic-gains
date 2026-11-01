function fig3_phase_paper(tags, SavePath)
% FIG3_PHASE_PAPER  论文 Fig.3 候选：把"分期"画在 (B,G) 参数平面上。
%
%   fig3_phase_paper()                        % 默认读 A = 3.25 / 4.00 / 4.60 三个平面
%   fig3_phase_paper(tags, SavePath)
%
% 输入
%   tags       scan_bg 输出后缀，如 {'A3p25','A4p00','A4p60'}（按列顺序）
%   SavePath   导出路径（矢量 PDF）；省略则只显示
%
% 输出
%   无返回值；有 SavePath 时导出 PDF
%
% 做什么
%   两行三列：
%     (a)(b)(c) 三个 A 值下的"主导期"分区图 —— 看区域位置怎么随 A 移动；
%     (d) 快活动占比 relFast 的连续场 + 阈值等值线（A 取最后一个 tag）；
%     (e) 主频 f_peak —— 与 (d) 互为独立证据（一个是能量比，一个是谱峰位置）；
%     (f) 四个区的代表波形各 5 s —— 让读者看到分区对应的实际 EEG 形态。
%
% 调用关系
%   上游：run_fig3_phase.m
%   下游：tests/results/scan_bg_<tag>.mat（由 tests/scan_bg.m 生成）、simulate_network
%
% 备注
%   分区图每个格点是该格点 50 s 内 21 个 10 s 窗的多数类；格点颜色只表达"主导期"，
%   格点内部的抖动（其余类占比）在 scan_bg 的表里，图上不表达。
%   (f) 的四条波形各自按自身峰值归一化（幅度差 20 倍，共用刻度就看不见小的），
%   实际幅度标在右侧。
if nargin < 1 || isempty(tags)
    tags = {'A3p25','A4p00','A4p60'};
end
if nargin < 2, SavePath = []; end

here = fileparts(mfilename('fullpath'));
resdir = fullfile(fileparts(here), 'tests', 'results');

%% ===== 全局绘图样式（与 fig1/fig3_paper 一致）=====
set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultLineLineWidth', 0.9);

% 类别配色：背景浅灰；其余按"发作间期→发作"的暖→冷顺序
cmap = [0.86 0.86 0.86;   % intICTAL
        0.96 0.78 0.45;   % preONSET
        0.85 0.45 0.75;   % ONSET
        0.25 0.50 0.80;   % ICTAL
        0.50 0.75 0.50;   % postICTAL
        0.60 0.60 0.60];  % uncertain

f = figure('Visible','off', 'Color','w', 'Units','inches', ...
    'Position',[1 1 6.3 3.9], 'PaperUnits','inches', ...
    'PaperSize',[6.3 3.9], 'PaperPosition',[0 0 6.3 3.9], ...
    'PaperPositionMode','manual');
tl = tiledlayout(2, 3, 'Padding','compact', 'TileSpacing','compact');

S = cell(1, numel(tags));
for k = 1:numel(tags)
    S{k} = load(fullfile(resdir, ['scan_bg_' tags{k} '.mat']));
end
names = S{1}.names;

%% ===== 上排：三个 A 下的分区图 =====
for k = 1:numel(tags)
    ax = nexttile;
    codes = codes_of(S{k}, names);
    imagesc(ax, S{k}.Gv, S{k}.Bv, codes);
    axis(ax, 'xy'); colormap(ax, cmap); clim(ax, [0.5 numel(names)+0.5]);
    set(ax, 'XTick', 0:10:50, 'YTick', 5:10:55, 'TickDir','out', 'Box','on', 'Layer','top');
    xlabel(ax, 'G_f (mV)'); ylabel(ax, 'G_s (mV)');
    title(ax, sprintf('(%s) G_e = %.2f mV', char('a'+k-1), S{k}.A0), 'FontWeight','normal');
    if k == numel(tags)
        cb = colorbar(ax); cb.Ticks = 1:numel(names); cb.TickLabels = names;
        cb.FontSize = 6.5; cb.Box = 'off';
    end
end

%% ===== 下排左二：连续场 =====
Sl = S{end};
ax = nexttile;
imagesc(ax, Sl.Gv, Sl.Bv, Sl.relfast_g); axis(ax,'xy');
colormap(ax, flipud(gray)); clim(ax, [0 0.9]);
cb = colorbar(ax); cb.FontSize = 6.5; cb.Box = 'off';
set(ax, 'XTick', 0:10:50, 'YTick', 5:10:55, 'TickDir','out', 'Box','on', 'Layer','top');
xlabel(ax, 'G_f (mV)'); ylabel(ax, 'G_s (mV)');
title(ax, '(d) fast-activity fraction, 15-40 Hz', 'FontWeight','normal');
hold(ax, 'on');
% 分类阈值等值线：读者能直接看到"颜色到哪一格会换标签"
contour(ax, Sl.Gv, Sl.Bv, Sl.relfast_g, [Sl.p.thr_fast Sl.p.thr_fast], ...
    'LineColor', [0.85 0.15 0.15], 'LineWidth', 1.3);
hold(ax, 'off');

ax = nexttile;
imagesc(ax, Sl.Gv, Sl.Bv, Sl.fpeak_g); axis(ax,'xy');
colormap(ax, parula); clim(ax, [0 30]);
cb = colorbar(ax); cb.FontSize = 6.5; cb.Box = 'off';
set(ax, 'XTick', 0:10:50, 'YTick', 5:10:55, 'TickDir','out', 'Box','on', 'Layer','top');
xlabel(ax, 'G_f (mV)'); ylabel(ax, 'G_s (mV)');
title(ax, '(e) peak frequency (Hz)', 'FontWeight','normal');

%% ===== 下排右：四个区的代表波形（三个小轴叠起来）=====
ax = nexttile;
ax.Visible = 'off';
pos = ax.Position;
which = [ 0 30;   % ONSET
         25 10;   % ICTAL
         45 20;   % intICTAL
         35 10];  % preONSET
lab_txt = {'ONSET','ICTAL','intICTAL','preONSET'};
Fs = 256; T = 5; n = size(which,1);
hgap = 0.085;
h = (pos(4) - (n-1)*hgap) / n;
for k = 1:n
    a = axes('Position', [pos(1), pos(2) + (n-k)*(h+hgap), pos(3), h]); %#ok<LAXES>
    cfg = default_config();
    cfg.node_count = 1; cfg.model_dynamic = 0; cfg.esf_enable = false;
    cfg.resample_dy = false; cfg.simulation_time = 40;
    cfg.A = Sl.A0; cfg.B = which(k,1); cfg.G = which(k,2);
    cfg.seed = 5150 + k;
    out = simulate_network(cfg);
    % 选段：默认取"峰均比最接近中位数"的 5 s 段（有代表性）；
    % 棘波区（第 4 行）取峰均比最大的段，否则很可能正巧落在棘波之间的平静期。
    seg = 5; k0 = round(28*Fs);
    rrs = []; idx = [];
    for s0 = round(22*Fs) : round(2*Fs) : (numel(out) - seg*Fs - 1)
        xs = out(s0+1 : s0+seg*Fs);
        rrs(end+1) = max(abs(xs)) / max(sqrt(mean(xs.^2)), eps); %#ok<SAGROW>
        idx(end+1) = s0;                                          %#ok<SAGROW>
    end
    if k == 4
        [~, ii] = max(rrs);
    else
        [~, ii] = min(abs(rrs - median(rrs)));
    end
    k0 = idx(ii);
    x = out(k0+1 : k0+T*Fs);
    amp = prctile(abs(x), 99);      % 用 p99 归一化：个别孤立大棘波不该把整条压平
    rms = sqrt(mean((x-mean(x)).^2));
    t = (0:numel(x)-1)/Fs;
    plot(a, t, x/amp, 'k');
    ylim(a, [-1.35 1.35]); xlim(a, [0 T]);
    set(a, 'YTick', [], 'XTick', [], 'Box','on', 'TickDir','out');
    text(a, -0.12, 0.5, lab_txt{k}, 'Units','normalized', ...
        'HorizontalAlignment','right', 'FontSize', 6.5);
    text(a, 1.01, 0.5, sprintf('rms %.1f mV', rms), 'Units','normalized', ...
        'FontSize', 6, 'Color', [0.35 0.35 0.35]);
    if k == n
        set(a, 'XTick', 0:1:T);
        xlabel(a, 'time (s)');
    end
end
% ax 已隐藏（只用来占位）。面板编号放成左侧小标签，避免与上排的 x 轴标签打架
annotation(f, 'textbox', [pos(1)-0.080, pos(2)+pos(4)+0.008, 0.05, 0.035], ...
    'String', '(f)', 'LineStyle','none', 'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', 'FontSize', 9, 'FontName','Arial');

if ~isempty(SavePath)
    print(f, SavePath, '-dpdf', '-r300');
end
close(f);
end

%% ================= 本地工具 =================
function codes = codes_of(S, names)
codes = zeros(numel(S.Bv), numel(S.Gv));
for c = 1:numel(names)
    codes(strcmp(S.lab_grid, names{c})) = c;
end
end
