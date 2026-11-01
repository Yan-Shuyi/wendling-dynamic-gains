function scan_bg(Bv, Gv, A0, sfx)
% SCAN_BG  (B,G) 平面上的静态增益参数扫描：每个格点判"这是哪个期"。
%
%   用途
%     回答导师的第二问（fig3 用时间分期是不是没意义）：把"分期"从时间轴搬到
%     参数轴上——如果同一套 v2 判据在 (B,G) 平面上划出成片的区域，那"期"就是
%     参数区的名字，而不是参数漂移过程中的时间产物；如果扫不出成片区域、只有
%     零星散点，那说明判据其实是在给噪声贴标签，fig3 就得改写法。
%
%   做法
%     1. 固定 A = 3.25（兴奋性增益），只扫慢抑制 B 与快抑制 G；
%     2. 静态增益（model_dynamic = 0）：每个格点是一个**定常**动力学系统，
%        没有 ESF 触发、没有增益漂移，所以"这个格点是什么期"是well-defined 的；
%     3. 每个格点跑 70 s，丢掉前 20 s 暂态，剩下 50 s 用 10 s 窗、2 s 步长切，
%        逐窗跑 classify_v2，取出现最多的类为该格点的类，并记录各类占比；
%     4. 幅度判据的参考段用**全平面的中位数**（不逐点自比，否则"比自己开头大
%        2.5 倍"这种循环定义会让每个格点都容易判成 ICTAL）。
%
%   输出
%     results/scan_bg_<A>.txt   逐格点表格（<A> 形如 A3p25）
%     results/scan_bg_<A>.mat    网格矩阵（画图/复查用）
%     figures/fig3_phase_BG_<A>.png  三张平面图：分类、快活动占比、主频
%
%   调用
%     addpath('C:/Users/26442/WorkBuddy/Claw/code'); scan_bg
%     scan_bg(5:10:55, 0:10:50)     % 也可指定网格（冒烟测试用）
%     scan_bg(0:2.5:15, 20:2.5:50, 4.60, 'fine')   % 细网格，输出名另存不覆盖粗网格
%
%   注意
%     静态模式下 System 的 dy_size = 0，ESF 跳变会撞维度错误——本脚本不开 ESF，
%     但 wendling.ESF 已按模式分流修好（见那里注释）。

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
addpath(fileparts(here));
[~, figdir] = setup_paths();

%% ---------------- 网格与仿真设置 ----------------
if nargin < 2 || isempty(Bv), Bv = 5:5:55; end   % 慢抑制增益 (mV)
if nargin < 1 || isempty(Gv), Gv = 0:5:50; end   % 快抑制增益 (mV)
if nargin < 3 || isempty(A0), A0 = 3.25; end     % 兴奋性增益 (mV)
if nargin < 4 || isempty(sfx), sfx = ''; end
tag = [strrep(sprintf('A%.2f', A0), '.', 'p') sfx];  % 输出后缀：A 值 + 可选备注（如 'fine'）
SIM_T   = 70;           % 单次仿真时长 (s)
T_TRANS = 20;           % 丢弃的暂态 (s)
p = classify_params_v2();
WIN = p.win_sec; STRIDE = 2;      % 10 s 窗、2 s 步长（比 1 s 省一半时间）

nB = numel(Bv); nG = numel(Gv);
EEG = cell(nB, nG);
lab_grid = cell(nB, nG);
rms_g = nan(nB, nG);  relfast_g = nan(nB, nG); relslow_g = nan(nB, nG);
fpeak_g = nan(nB, nG); spk_g = nan(nB, nG);
frac = zeros(nB, nG, 6);         % 六个类别的占比
names = p.names;

fprintf('\n===== (B,G) 静态扫描：%d × %d = %d 个格点，每点 %.0f s =====\n', ...
    nB, nG, nB*nG, SIM_T);
tic;
for i = 1:nB
    for j = 1:nG
        cfg = default_config();
        cfg.node_count    = 1;
        cfg.model_dynamic = 0;          % 静态增益：格点 = 一个定常系统
        cfg.esf_enable    = false;      % 不触发，纯看该增益下的自持行为
        cfg.resample_dy   = false;
        cfg.simulation_time = SIM_T;
        cfg.A = A0;   cfg.B = Bv(i); cfg.G = Gv(j);
        cfg.seed = 1000 + 10*i + j;     % 逐点固定种子，可复现
        try
            out = simulate_network(cfg);
            EEG{i,j} = out(:)';
        catch err
            fprintf(2, '  B=%g G=%g 仿真失败：%s\n', Bv(i), Gv(j), err.message);
            EEG{i,j} = [];
        end
    end
    fprintf('  B = %2.0f 完成（%.1f s）\n', Bv(i), toc);
end

%% ---------------- 全平面参考段标定 ----------------
% 参考 = 全平面 RMS / 斜率 p99 的中位数。理由：平面里绝大部分格点是"背景"，
% 中位数就是背景本身；若用某个格点自比，ICTAL 的幅度门槛就退化成
% "比自己开头大 2.5 倍"，静态系统里没有"开头"可言。
Fs = 256;
allrms = []; allslope = [];
for i = 1:nB
    for j = 1:nG
        x = EEG{i,j};
        if isempty(x), continue; end
        x2 = x(round(T_TRANS*Fs)+1 : end);
        allrms(end+1) = sqrt(mean(highpass(x2, p.hp_cut, Fs).^2));
        d = slope_series(x2, Fs, p.spk_lp);
        allslope(end+1) = prctile(d, 99);
    end
end
ref.rms = median(allrms);
ref.slope_p99 = median(allslope);
ref.range = ref.rms;      % 本脚本不用
ref.psd = []; ref.f = [];
fprintf('参考段（全平面中位数）：rms = %.3f mV，slope_p99 = %.1f mV/s\n', ...
    ref.rms, ref.slope_p99);

%% ---------------- 逐格点分期 ----------------
NW = floor((SIM_T - T_TRANS - WIN)/STRIDE) + 1;
for i = 1:nB
    for j = 1:nG
        x = EEG{i,j};
        if isempty(x), lab_grid{i,j} = 'failed'; continue; end
        rf = []; rs = []; fp = []; sp = []; labs = cell(1,NW);
        for k = 1:NW
            a = round((T_TRANS + (k-1)*STRIDE)*Fs) + 1;
            b = a + WIN*Fs - 1;
            if b > numel(x), break; end
            [lb, ~, f] = classify_v2(x(a:b), Fs, p, ref);
            labs{k} = lb;
            rf(end+1) = f.rel_fast; rs(end+1) = f.rel_slow;
            fp(end+1) = f.f_peak;   sp(end+1) = f.spk_rate;
        end
        labs = labs(~cellfun(@isempty, labs));
        for c = 1:numel(names)
            frac(i,j,c) = mean(strcmp(labs, names{c}));
        end
        % 该格点的类 = 出现最多的类（平票时按判据严格度取靠前的）
        cnt = squeeze(frac(i,j,:));
        [~, ic] = max(cnt);
        lab_grid{i,j} = names{ic};
        rms_g(i,j)     = sqrt(mean(highpass(x, p.hp_cut, Fs).^2));
        relfast_g(i,j) = median(rf);
        relslow_g(i,j) = median(rs);
        fpeak_g(i,j)   = median(fp);
        spk_g(i,j)     = max(sp);
    end
end

%% ---------------- 写表 ----------------
txtfile = fullfile(outdir, ['scan_bg_' tag '.txt']);
fid = fopen(txtfile, 'w');
fprintf(fid, '(B,G) 静态扫描：A = %.2f mV，静态增益，70 s，丢弃前 20 s；窗 10 s / 步 2 s\n', A0);
fprintf(fid, '参考段 = 全平面中位数：rms %.3f mV，slope_p99 %.1f mV/s\n', ref.rms, ref.slope_p99);
fprintf(fid, '判据：v2（Wendling 2005 四期，见 features/classify_params_v2.m）\n\n');
fprintf(fid, '%5s %5s %9s %9s %9s %8s %8s  %-11s %s\n', 'B', 'G', 'rms', 'relFast', 'relSlow', ...
    'fPeak', 'spkMax', 'label', '各期占比(i/p/O/I/u)');
for i = 1:nB
    for j = 1:nG
        c = squeeze(frac(i,j,:));
        fprintf(fid, '%5.0f %5.0f %9.3f %9.3f %9.3f %8.1f %8.2f  %-11s %.2f/%.2f/%.2f/%.2f/%.2f\n', ...
            Bv(i), Gv(j), rms_g(i,j), relfast_g(i,j), relslow_g(i,j), fpeak_g(i,j), spk_g(i,j), ...
            lab_grid{i,j}, c(1), c(2), c(3), c(4), c(6));
    end
end
fclose(fid);

%% ---------------- 汇总 ----------------
fprintf('\n各类格点数：\n');
for c = 1:numel(names)
    n = sum(sum(strcmp(lab_grid, names{c})));
    fprintf('  %-11s %3d / %d\n', names{c}, n, nB*nG);
end
fprintf('最大 relFast = %.3f（阈值 %.2f）出现在 B=%g, G=%g\n', ...
    max(relfast_g(:)), p.thr_fast, ...
    Bv(rem(find(relfast_g(:)==max(relfast_g(:)),1)-1, nB)+1), ...
    Gv(floor((find(relfast_g(:)==max(relfast_g(:)),1)-1)/nB)+1));

save(fullfile(outdir, ['scan_bg_' tag '.mat']), 'Bv', 'Gv', 'A0', 'lab_grid', 'frac', 'rms_g', ...
    'relfast_g', 'relslow_g', 'fpeak_g', 'spk_g', 'ref', 'p', 'names');

%% ---------------- 画平面图 ----------------
codes = zeros(nB, nG);
for c = 1:numel(names)
    codes(strcmp(lab_grid, names{c})) = c;
end
cmap = [0.85 0.85 0.85;   % intICTAL
        0.95 0.75 0.35;   % preONSET
        0.85 0.40 0.85;   % ONSET
        0.30 0.55 0.85;   % ICTAL
        0.45 0.75 0.45;   % postICTAL
        0.55 0.55 0.55];  % uncertain

f = figure('Visible','off','Position',[100 100 1250 380]);
subplot(1,3,1);
imagesc(Gv, Bv, codes); axis xy; colormap(gca, cmap); clim([0.5 numel(names)+0.5]);
cb = colorbar; cb.Ticks = 1:numel(names); cb.TickLabels = names;
xlabel('G（快抑制增益, mV）'); ylabel('B（慢抑制增益, mV）');
title('(a) 各格点的主导期');

subplot(1,3,2);
imagesc(Gv, Bv, relfast_g); axis xy; colorbar;
xlabel('G (mV)'); ylabel('B (mV)');
title(sprintf('(b) 快活动占比 P_{15-40}/P_{1-90}（阈值 %.2f）', p.thr_fast));

subplot(1,3,3);
imagesc(Gv, Bv, fpeak_g); axis xy; colorbar;
xlabel('G (mV)'); ylabel('B (mV)'); title('(c) 主频 (Hz)');
pngfile = fullfile(figdir, ['fig3_phase_BG_' tag '.png']);
print(f, pngfile, '-dpng', '-r150');
close(f);
fprintf('表：%s\n图：%s\n', txtfile, pngfile);
end
