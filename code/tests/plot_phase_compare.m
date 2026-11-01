function plot_phase_compare(tags)
% PLOT_PHASE_COMPARE  把多个 A 值下的 (B,G) 平面拼成一张对比图（论文 fig3 候选）。
%
%   plot_phase_compare()                 % 默认 A = 3.25 / 4.00 / 4.60
%   plot_phase_compare({'A3p25','A4p60'})
%
% 输入
%   tags   scan_bg 输出的文件后缀列表；每个对应 results/scan_bg_<tag>.mat
%
% 输出
%   figures/fig3_phase_compare.png   行 = A，列 = 分类图 / 快活动占比 / 主频
%
% 说明
%   数据全部来自 scan_bg（静态增益扫描），本函数只负责排版，不重新仿真。
%   分类色块之间用灰阶区分，跨 A 对比时看"色块的形状与位置怎么移动"。
if nargin < 1 || isempty(tags)
    tags = {'A3p25','A4p00','A4p60'};
end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
[~, figdir] = setup_paths();
resdir = fullfile(here, 'results');

nR = numel(tags);
cmap = [0.85 0.85 0.85;   % intICTAL
        0.95 0.75 0.35;   % preONSET
        0.85 0.40 0.85;   % ONSET
        0.30 0.55 0.85;   % ICTAL
        0.45 0.75 0.45;   % postICTAL
        0.55 0.55 0.55];  % uncertain

f = figure('Visible','off','Position',[80 80 1250 300*nR]);
for r = 1:nR
    S = load(fullfile(resdir, ['scan_bg_' tags{r} '.mat']));
    names = S.names; Bv = S.Bv; Gv = S.Gv;
    codes = zeros(numel(Bv), numel(Gv));
    for c = 1:numel(names)
        codes(strcmp(S.lab_grid, names{c})) = c;
    end

    subplot(nR, 3, 3*r-2);
    imagesc(Gv, Bv, codes); axis xy; colormap(gca, cmap);
    clim([0.5 numel(names)+0.5]);
    cb = colorbar; cb.Ticks = 1:numel(names); cb.TickLabels = names;
    ylabel('B（慢抑制, mV）');
    if r == nR, xlabel('G（快抑制, mV）'); end
    title(sprintf('A = %.2f：各格点主导期', S.A0));

    subplot(nR, 3, 3*r-1);
    imagesc(Gv, Bv, S.relfast_g); axis xy; colorbar; clim([0 0.9]);
    ylabel('B (mV)');
    if r == nR, xlabel('G (mV)'); end
    title(sprintf('A = %.2f：P_{15-40}/P_{1-90}', S.A0));

    subplot(nR, 3, 3*r);
    imagesc(Gv, Bv, S.fpeak_g); axis xy; colorbar; clim([0 30]);
    ylabel('B (mV)');
    if r == nR, xlabel('G (mV)'); end
    title(sprintf('A = %.2f：主频 (Hz)', S.A0));
end
out = fullfile(figdir, 'fig3_phase_compare.png');
print(f, out, '-dpng', '-r150');
close(f);
fprintf('图：%s\n', out);
end
