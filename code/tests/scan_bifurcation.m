function scan_bifurcation(stage)
% SCAN_BIFURCATION  一维参数扫描：把"期"读成参数的分岔结构，而不是我们贴的标签。
%
%   scan_bifurcation()        % 两段都跑
%   scan_bifurcation('S1')    % 只跑 S1（Gs 扫描，A = 4.60，两条 Gf）
%   scan_bifurcation('S2')    % 只跑 S2（Ge 扫描，两组抑制参数）
%
% 两段扫描
%   [S1] 慢抑制 Gs = 0:1:50，A = 4.60（发作区存在的兴奋性水平），Gf = 10 与 30。
%        Gf = 10 时沿 Gs 下降应当从发作区退出到背景；Gf = 30 时应当从发作区进入快活动区。
%        这两条是同一个"崩塌"事件在两种 Gf 下的两种落点——Fig.7 的主面板。
%   [S2] 兴奋性 Ge = 2.5:0.1:6.0，取 (Gs,Gf) = (25,10) 与 (10,10)。
%        看"慢节律在 Ge 超过某个值之后才出现"。
%
% 做法
%   静态增益（无 ESF），每个参数点跑 70 s、丢前 20 s，10 s 窗 / 2 s 步长逐窗记
%   主频、rms、快/慢带占比与 v2 标签。每个窗是一个点，所以同一参数会竖着排一列点。
%
% 输出
%   results/bifur_<stage>.txt / .mat
%
% 耗时：S1 约 6 min，S2 约 6 min。

if nargin < 1 || isempty(stage), stage = 'both'; end
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();

SIM_T = 70; T_TRANS = 20;
p = classify_params_v2();
WIN = p.win_sec; STRIDE = 2;
Fs = 256;

%% ---------------- S1：Gs 扫描（A = 4.60） ----------------
if any(strcmpi(stage, {'both','S1'}))
    S.x = 0:1:50;
    S.name = 'Gs';
    S.A = 4.60;
    S.sets = [10; 30];                 % 两种 Gf
    S.lbl = {'Gf=10', 'Gf=30'};
    S.note = 'A = 4.60 mV；Gf = 10 与 30';

    fprintf('\n[S1] Gs 扫描：%d 点 × %d 条 = %d 次\n', numel(S.x), size(S.sets,1), ...
        numel(S.x)*size(S.sets,1));
    tic;
    C = run_sweep(S, SIM_T, 5000);
    fprintf('  仿真 %.0f s\n', toc);
    F = score(C, T_TRANS, WIN, STRIDE, Fs, p);
    save(fullfile(outdir, 'bifur_S1.mat'), 'S', 'F');
    write_txt(fullfile(outdir, 'bifur_S1.txt'), S, F, 'S1', SIM_T, T_TRANS);
end

%% ---------------- S2：Ge 扫描 ----------------
if any(strcmpi(stage, {'both','S2'}))
    S.x = 2.5:0.1:6.0;
    S.name = 'Ge';
    S.A = [];                          % Ge 就是自变量
    S.sets = [25 10; 10 10];           % [Gs Gf]
    S.lbl = {'Gs=25,Gf=10', 'Gs=10,Gf=10'};
    S.note = '(Gs,Gf) = (25,10) 与 (10,10)';

    fprintf('\n[S2] Ge 扫描：%d 点 × %d 条 = %d 次\n', numel(S.x), size(S.sets,1), ...
        numel(S.x)*size(S.sets,1));
    tic;
    C = run_sweep(S, SIM_T, 7000);
    fprintf('  仿真 %.0f s\n', toc);
    F = score(C, T_TRANS, WIN, STRIDE, Fs, p);
    save(fullfile(outdir, 'bifur_S2.mat'), 'S', 'F');
    write_txt(fullfile(outdir, 'bifur_S2.txt'), S, F, 'S2', SIM_T, T_TRANS);
end

% 参考段需要两段合并才能算，单独跑一段时用同一条轨迹内部的中位数
fprintf('\n完成。\n');
end

% ----------------------------------------------------------------------
function C = run_sweep(S, SIM_T, seed0)
C = cell(numel(S.x), size(S.sets, 1));
for j = 1:size(S.sets, 1)
    for i = 1:numel(S.x)
        cfg = default_config();
        cfg.node_count = 1; cfg.model_dynamic = 0; cfg.esf_enable = false;
        cfg.resample_dy = false; cfg.simulation_time = SIM_T;
        if isempty(S.A)                    % Ge 扫描
            cfg.A = S.x(i);  cfg.B = S.sets(j,1);  cfg.G = S.sets(j,2);
        else                               % Gs 扫描
            cfg.A = S.A;     cfg.B = S.x(i);       cfg.G = S.sets(j,1);
        end
        cfg.seed = seed0 + 100*j + i;
        try
            out = simulate_network(cfg);
            C{i, j} = out(:)';
        catch err
            fprintf(2, '  失败 %s=%g：%s\n', S.name, S.x(i), err.message);
        end
    end
    fprintf('  %s 完成\n', S.lbl{j});
end
end

% ----------------------------------------------------------------------
function F = score(C, T_TRANS, WIN, STRIDE, Fs, p)
% 参考段：本段扫描全部格点的中位数（静态系统没有"自己的开头"，自比会让幅度门槛循环）
allrms = []; allslope = [];
for k = 1:numel(C)
    x = C{k};
    if isempty(x), continue; end
    x2 = x(round(T_TRANS*Fs)+1 : end);
    allrms(end+1)   = sqrt(mean(highpass(x2, p.hp_cut, Fs).^2));      %#ok<SAGROW>
    allslope(end+1) = prctile(slope_series(x2, Fs, p.spk_lp), 99);    %#ok<SAGROW>
end
ref.rms = median(allrms);  ref.slope_p99 = median(allslope);
ref.range = ref.rms;  ref.psd = [];  ref.f = [];
fprintf('  参考段：rms %.3f mV，slope_p99 %.1f mV/s\n', ref.rms, ref.slope_p99);

[nx, ns] = size(C);
NW = floor((numel(C{1})/Fs - T_TRANS - WIN)/STRIDE) + 1;
F.f = nan(nx,ns,NW);  F.rms = nan(nx,ns,NW);
F.rf = nan(nx,ns,NW); F.rs = nan(nx,ns,NW); F.spk = nan(nx,ns,NW);
F.lab = cell(nx,ns,NW);
for j = 1:ns
    for i = 1:nx
        x = C{i,j};
        if isempty(x), continue; end
        for w = 1:NW
            t0 = T_TRANS + (w-1)*STRIDE;
            k1 = round(t0*Fs) + 1; k2 = k1 + round(WIN*Fs) - 1;
            if k2 > numel(x), continue; end
            [lab, ~, f] = classify_v2(x(k1:k2), Fs, p, ref);
            F.lab{i,j,w} = lab;
            F.f(i,j,w) = f.f_peak;  F.rms(i,j,w) = f.rms;
            F.rf(i,j,w) = f.rel_fast;  F.rs(i,j,w) = f.rel_slow;
            F.spk(i,j,w) = f.spk_rate;
        end
    end
end
end

% ----------------------------------------------------------------------
function write_txt(fname, S, F, tag, SIM_T, T_TRANS)
fid = fopen(fname, 'w');
fprintf(fid, '%s：%s 扫描（%s），每点 %.0f s（丢前 %.0f s），10 s 窗 / 2 s 步长\n', ...
    tag, S.name, S.note, SIM_T, T_TRANS);
for j = 1:size(S.sets, 1)
    fprintf(fid, '--- %s ---\n%8s %8s %8s %8s %8s %8s %8s  %s\n', S.lbl{j}, S.name, ...
        'f_med', 'f_lo', 'f_hi', 'rms', 'relS', 'relF', 'dominant');
    for i = 1:numel(S.x)
        f  = squeeze(F.f(i,j,:));   rm = squeeze(F.rms(i,j,:));
        rs = squeeze(F.rs(i,j,:));  rf = squeeze(F.rf(i,j,:));
        lb = squeeze(F.lab(i,j,:));
        if all(isnan(f)), continue; end
        [u, ~, ic] = unique(lb);
        c = accumarray(ic, 1);  [~, im] = max(c);
        fprintf(fid, '%8.2f %8.1f %8.1f %8.1f %8.2f %8.2f %8.2f  %s\n', S.x(i), ...
            median(f,'omitnan'), min(f), max(f), median(rm,'omitnan'), ...
            median(rs,'omitnan'), median(rf,'omitnan'), u{im});
    end
end
fclose(fid);
fprintf('  写出 %s\n', fname);
end
