function event_stats(nodes, ncls, ges)
% EVENT_STATS  事件对齐的分期统计：把"四个固定时刻"换成"相对 onset 的四个相位"。
%
%   event_stats                       % 默认 60 例 × 3 个触发幅度
%   event_stats(60, 20, [4.00 4.25 4.50])
%
% 为什么重做
%   §4.3 原来在 10 / 30 / 50 / 70 s 四个"代表时刻"上取增益画箱线图。D14 已经证明这件事不稳：
%   触发幅度只要动几个百分点，onset 就整体平移，固定时刻取到的其实是不同相位，Gs 的 IQR
%   可以从 0.2 涨到 20 以上。根因是网络里没有统一的钟——每例的 onset 落在不同时刻。
%   这里对每个触发幅度各跑一批，两种取法都算：固定时刻 vs 按 onset 对齐，
%   再用"跨触发幅度的一致性"把差别量化出来（同一相位、不同幅度，中位数应当一致）。
%
% 同时回答另一件事
%   §4.2 写了"扰动前标签在 interictal 与发作之间交替、增益进入发作枝后转为 ictal"。
%   这句话当时是按 §36 的结论写的，没有实测。这里对若干例逐窗跑 classify_v2，
%   统计"扰动前 / onset 之后"两类时段的标签构成，用数据钉死或推翻它。
%
% 输出
%   results/event_stats.txt   表格（对齐 vs 固定、跨幅度一致性、标签构成）
%   results/event_stats.mat   原始数组（画图与复查用）
%
% 耗时：3 × 60 例 × 200 s 动态仿真约 9 min；标签统计 20 例 × 96 窗约 1 min。

if nargin < 1 || isempty(nodes), nodes = 60; end
if nargin < 2 || isempty(ncls),  ncls  = 20; end
if nargin < 3 || isempty(ges),   ges   = [4.00 4.25 4.50]; end

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();

% 四个相位的定义：发作间期用绝对时刻（扰动前，各幅度相同），其余相对 onset。
% 不用更早的相对相位（如 -30 s）：onset 早的例子那里还在代偿上升段，
% onset 晚的例子那里还没扰动，取到的本来就不是同一状态。
ph_kind = {'abs', 'rel', 'rel', 'rel'};
ph_val  = [9, -10, 0, 20];
name_al = {'interictal', 'pre-onset(-10)', 'onset(0)', 'ictal(+20)'};
t_fix   = [9 30 50 70];
name_fx = {'9s', '30s', '50s', '70s'};
nph     = numel(ph_val);
gn      = {'Ge', 'Gs', 'Gf'};
nA      = numel(ges);

% 存每批结果
T_all    = cell(1, nA);          % onset
al_all   = cell(nA, 3);          % {幅度}(gain) = 相位 × 例
fx_all   = cell(nA, 3);
traj     = cell(1, nA);          % 对齐后的 Gs 轨迹（供画图），nTraj × 161
tg       = -40:0.5:40;
cfg0     = default_config();
cfg0.noise_segments = [0.3 0.5 0.2];   % 历史设置：Fig.5 的存量结果基于三段噪声

for a = 1:nA
    cfg = cfg0;
    cfg.node_count = nodes;
    cfg.A          = 3.25;
    cfg.B          = 22;
    cfg.G          = 10;
    cfg.esf_enable = true;
    cfg.esf_time   = 10;
    cfg.esf_params = [ges(a) 22 10];
    cfg.esf_jitter = 0.05;
    cfg.resample_dy = true;

    fprintf('event_stats: 触发 Ge=%.2f，%d 例 × %.0f s …\n', ges(a), nodes, cfg.simulation_time);
    tic;
    [out_fs_a{a}, dy_fs, cfg] = simulate_network(cfg); %#ok<AGROW>
    fprintf('  %.1f s\n', toc);

    fs = cfg.dy_fs;  n = size(dy_fs, 2);  t = (0:n-1) / fs;
    Ge = dy_fs(1:3:end, :);  Gs = dy_fs(2:3:end, :);  Gf = dy_fs(3:3:end, :);
    G  = {Ge, Gs, Gf};

    % onset：Gs 首次跌破中点（静息 22、下限约 3.5 的中点）
    thr = 12.75;
    T = NaN(nodes, 1);
    for i = 1:nodes
        k = find(Gs(i, :) < thr, 1, 'first');
        if ~isempty(k), T(i) = t(k); end
    end
    ok = ~isnan(T);
    T_all{a} = T;
    fprintf('  onset 检出 %d/%d，中位 %.1f s（%.1f–%.1f）\n', sum(ok), nodes, ...
        median(T(ok)), min(T(ok)), max(T(ok)));

    for g = 1:3
        Aa = NaN(nph, nodes);
        Af = NaN(numel(t_fix), nodes);
        for q = 1:nph
            for i = find(ok)'
                if strcmp(ph_kind{q}, 'abs'), tq = ph_val(q); else, tq = T(i) + ph_val(q); end
                k = round(tq * fs) + 1;
                if k >= 1 && k <= n, Aa(q, i) = G{g}(i, k); end
            end
        end
        for q = 1:numel(t_fix)
            k = round(t_fix(q) * fs) + 1;
            Af(q, :) = G{g}(:, k)';
        end
        al_all{a, g} = Aa;
        fx_all{a, g} = Af;
    end

    % 供画图的对齐轨迹（Gs），落在公共时间栅格上
    S = NaN(nodes, numel(tg));
    for i = find(ok)'
        for q = 1:numel(tg)
            k = round((T(i) + tg(q)) * fs) + 1;
            if k >= 1 && k <= n, S(i, q) = Gs(i, k); end
        end
    end
    traj{a} = S;
end

%% ---------------- 表 1：每批的对齐 / 固定时刻 ----------------
fid = fopen(fullfile(outdir, 'event_stats.txt'), 'w');
fprintf(fid, '事件对齐 vs 固定时刻：每个触发幅度 %d 例，onset = Gs < 12.75 mV 首次跨越\n\n', nodes);

summary = struct();
for a = 1:nA
    fprintf(fid, '========== 触发 Ge = %.2f mV ==========\n', ges(a));
    T = T_all{a};  okT = ~isnan(T);
    fprintf(fid, 'onset 中位 %.1f s（%.1f–%.1f），检出 %d/%d\n', median(T(okT)), ...
        min(T(okT)), max(T(okT)), sum(okT), numel(T));
    for g = 1:3
        fprintf(fid, '--- %s ---\n', gn{g});
        [mi_a, sp_a, med_a] = dump_table(fid, al_all{a, g}, name_al, 'aligned');
        [mi_f, sp_f, med_f] = dump_table(fid, fx_all{a, g}, name_fx, 'fixed  ');
        summary.iqr_a(a, g) = mi_a;  summary.sep_a(a, g) = sp_a;  summary.med_a(a, g, :) = med_a;
        summary.iqr_f(a, g) = mi_f;  summary.sep_f(a, g) = sp_f;  summary.med_f(a, g, :) = med_f;
    end
    fprintf(fid, '\n');
end

%% ---------------- 表 2：跨触发幅度的一致性（这一张才是判决性的） ----------------
% 同名的"相位"在三种触发幅度下应当给出同一组增益统计量——如果这个"相位"真是同一相位。
% 对齐后成立；固定时刻不成立，且偏差随幅度增大而增大。
fprintf(fid, '========== 跨触发幅度一致性 ==========\n');
fprintf(fid, '同一标签的相位在 %d 个触发幅度之间的中位数极差（mV）\n', nA);
fprintf(fid, '%-5s %-16s %10s %10s\n', 'gain', 'phase', 'aligned', 'fixed');
for g = 1:3
    for q = 1:nph
        v_a = squeeze(summary.med_a(:, g, q));
        fprintf(fid, '%-5s %-16s %10.2f', gn{g}, name_al{q}, max(v_a) - min(v_a));
        if q <= numel(name_fx)
            v_f = squeeze(summary.med_f(:, g, q));
            fprintf(fid, ' %10.2f', max(v_f) - min(v_f));
        else
            fprintf(fid, ' %10s', '—');
        end
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n（fixed 列按序号对照：interictal↔9s、pre-onset↔30s、onset↔50s、ictal↔70s；\n');
fprintf(fid, ' 两者不是同一相位，这一列要读的是"它以为自己在取同一相位"——所以偏差就是错的相位）\n\n');

fprintf(fid, '=== 汇总：每批 max IQR (mV) ===\n');
fprintf(fid, '%-5s %10s %10s %10s\n', 'gain', 'aligned', 'fixed', 'ratio');
for g = 1:3
    fprintf(fid, '%-5s %10.2f %10.2f %10.1fx\n', gn{g}, max(summary.iqr_a(:, g)), ...
        max(summary.iqr_f(:, g)), max(summary.iqr_f(:, g)) / max(max(summary.iqr_a(:, g)), eps));
end

%% ---------------- 表 3：标签构成 ----------------
% 只对基准幅度做（它与 §4.2 的用例一致）
a0 = find(abs(ges - 4.25) < 1e-9, 1);
if isempty(a0), a0 = 1; end
cfg = cfg0; cfg.esf_time = 10;
out_fs = out_fs_a{a0};  T = T_all{a0};
fprintf(fid, '\n========== 分类器标签构成（触发 Ge=%.2f，%d 例，10 s 窗 / 2 s 步）==========\n', ...
    ges(a0), ncls);
p = classify_params_v2();
labels = {'intICTAL', 'preONSET', 'ONSET', 'ICTAL', 'uncertain'};

ref_x = [];
for i = 1:ncls
    ref_x = [ref_x, out_fs(i, 1:round(40 * cfg.eeg_fs))]; %#ok<AGROW>
end
ref = make_reference(ref_x, cfg.eeg_fs, p);

edges = [-40 -20 -10 0 10 20 40];
nbin  = numel(edges) - 1;
cnt   = zeros(nbin, numel(labels));
sub   = zeros(ncls, numel(labels), 2);
win = 10; step = 2;

for i = 1:ncls
    for t0 = 0:step:(cfg.simulation_time - win)
        k1 = round(t0 * cfg.eeg_fs) + 1;
        k2 = k1 + round(win * cfg.eeg_fs) - 1;
        lab = classify_v2(out_fs(i, k1:k2), cfg.eeg_fs, p, ref);
        li = find(strcmp(labels, lab), 1);
        if isempty(li), li = numel(labels); end
        tc = t0 + win/2;  rel = tc - T(i);
        j = find(rel >= edges(1:end-1) & rel < edges(2:end), 1, 'first');
        if ~isempty(j), cnt(j, li) = cnt(j, li) + 1; end
        if tc < cfg.esf_time
            sub(i, li, 1) = sub(i, li, 1) + 1;
        elseif rel > 5
            sub(i, li, 2) = sub(i, li, 2) + 1;
        end
    end
end

fprintf(fid, '相对相位 bin（相对 onset，s）：\n%-14s', 'bin');
for q = 1:numel(labels), fprintf(fid, ' %10s', labels{q}); end
fprintf(fid, '\n');
for j = 1:nbin
    fprintf(fid, '[%3d,%3d)    ', edges(j), edges(j+1));
    s = sum(cnt(j, :));
    for q = 1:numel(labels), fprintf(fid, ' %9.1f%%', 100 * cnt(j, q) / max(s, 1)); end
    fprintf(fid, '   (n=%d)\n', s);
end
fprintf(fid, '\n粗时段构成：\n');
for w = 1:2
    s = sum(sub(:, :, w), 1);  tot = sum(s);
    if w == 1, fprintf(fid, '  扰动前（窗中心 < %.0f s）  ', cfg.esf_time);
    else,      fprintf(fid, '  onset+5 s 之后          '); end
    for q = 1:numel(labels)
        fprintf(fid, '%s %.0f%%  ', labels{q}, 100 * s(q) / max(tot, 1));
    end
    fprintf(fid, '(n=%d)\n', tot);
end
fclose(fid);

save(fullfile(outdir, 'event_stats.mat'), 'ges', 'ph_kind', 'ph_val', 't_fix', 'name_al', 'name_fx', ...
    'T_all', 'al_all', 'fx_all', 'summary', 'cnt', 'edges', 'labels', 'sub', 'tg', 'traj');

fprintf('已写出：%s\n', fullfile(outdir, 'event_stats.txt'));
end

% ----------------------------------------------------------------------
function [max_iqr, sep, med] = dump_table(fid, A, names, tag)
med = median(A, 2, 'omitnan');
iqr = prctile(A, 75, 2) - prctile(A, 25, 2);
fprintf(fid, '  %s  %-16s %8s %8s %8s %8s\n', tag, 'phase', 'median', 'IQR', 'p5', 'p95');
for q = 1:numel(names)
    v = A(q, ~isnan(A(q, :)));
    fprintf(fid, '  %s  %-16s %8.2f %8.2f %8.2f %8.2f\n', tag, names{q}, ...
        med(q), iqr(q), prctile(v, 5), prctile(v, 95));
end
max_iqr = max(iqr);
d = abs(diff(med));
sep = min(d) / max(mean(iqr), eps);
fprintf(fid, '  %s  max IQR = %.2f，最小相邻相位中位数差 = %.2f，分离度 = %.2f\n', ...
    tag, max_iqr, min(d), sep);
end
