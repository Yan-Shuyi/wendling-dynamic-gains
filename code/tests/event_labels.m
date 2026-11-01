function event_labels(nodes)
% EVENT_LABELS  动态增益轨迹上的逐窗标签：参考段取"扰动前 10 s"，而不是整段开头。
%
%   event_labels()        % 默认 30 例
%   event_labels(30)
%
% 为什么单独做
%   event_stats 里为了统一口径，参考段用了 0–40 s 的拼接。对带增益动力学的轨迹这不合适：
%   扰动后增益爬升、EEG 幅度在 15–45 s 之间本身就抬高了一个量级，把这段算进参考，
%   等于把"发作"当成了基线，于是所有窗都判成 intICTAL——那是参考段的问题，不是模型的。
%   这里改成每例用自己的前 10 s（此时 Ge 尚未跳变、系统仍在静息）作参考。
%
% 输出
%   results/event_labels.txt   相对 onset 的标签构成（含三个粗时段）
%   results/event_labels.mat   逐窗标签（画图用）
%
% 耗时：30 例 × 200 s 约 2.5 min。

if nargin < 1 || isempty(nodes), nodes = 30; end
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();

cfg = default_config();
cfg.noise_segments = [0.3 0.5 0.2];   % 历史设置：Fig.5(d) 的存量结果基于三段噪声
cfg.node_count = nodes;
cfg.A = 3.25; cfg.B = 22; cfg.G = 10;
cfg.esf_enable = true; cfg.esf_time = 10;
cfg.esf_params = [4.25 22 10]; cfg.esf_jitter = 0.05;
cfg.resample_dy = true;

fprintf('event_labels: %d 例 × %.0f s …\n', cfg.node_count, cfg.simulation_time);
tic;
[out_fs, dy_fs, cfg] = simulate_network(cfg);
fprintf('  %.1f s\n', toc);

fs = cfg.dy_fs; n = size(dy_fs, 2); t = (0:n-1)/fs;
Gs = dy_fs(2:3:end, :);
thr = 12.75;
T = NaN(nodes, 1);
for i = 1:nodes
    k = find(Gs(i, :) < thr, 1, 'first');
    if ~isempty(k), T(i) = t(k); end
end

p = classify_params_v2();
labs = p.names;                       % 以 classify_params_v2 为准，避免名单漂移
nl = numel(labs);
win = 10; stride = 2;
Fs = cfg.eeg_fs;

edges = [-30 -10 0 10 20 40 80 160];
nbin = numel(edges) - 1;
cnt = zeros(nbin, nl);
coarse_name = {'pre-perturbation (<10 s)', 'build-up (10 s -> onset)', ...
               'after collapse (onset+5 -> onset+60)', 'late (onset+80 -> end)'};
csum = zeros(1, nl); cn = 0;
allab = cell(nodes, 1); atc = cell(nodes, 1);

for i = 1:nodes
    ref = make_reference(out_fs(i, 1:round(10*Fs)), Fs, p);
    labs_i = {}; tc_i = [];
    for t0 = 0:stride:(cfg.simulation_time - win)
        k1 = round(t0*Fs) + 1; k2 = k1 + round(win*Fs) - 1;
        lab = classify_v2(out_fs(i, k1:k2), Fs, p, ref);
        tc = t0 + win/2;
        labs_i{end+1} = lab; tc_i(end+1) = tc; %#ok<SAGROW>
        rel = tc - T(i);
        j = find(rel >= edges(1:end-1) & rel < edges(2:end), 1, 'first');
        li = find(strcmp(labs, lab), 1); if isempty(li), li = nl; end
        if ~isempty(j), cnt(j, li) = cnt(j, li) + 1; end
        if tc < cfg.esf_time || rel > 80
            csum(li) = csum(li) + 1;  cn = cn + 1;
        end
    end
    allab{i} = labs_i; atc{i} = tc_i;
end

fid = fopen(fullfile(outdir, 'event_labels.txt'), 'w');
fprintf(fid, '逐窗标签：%d 例，参考段 = 各例扰动前 10 s，10 s 窗 / 2 s 步\n', nodes);
fprintf(fid, 'onset 中位 %.1f s（%.1f–%.1f），检出 %d/%d\n\n', median(T(~isnan(T))), ...
    min(T(~isnan(T))), max(T(~isnan(T))), sum(~isnan(T)), nodes);
fprintf(fid, '%-14s', 'rel-onset bin');
for q = 1:nl, fprintf(fid, ' %9s', labs{q}); end
fprintf(fid, '\n');
for j = 1:nbin
    s = sum(cnt(j, :));
    fprintf(fid, '[%4d,%4d)   ', edges(j), edges(j+1));
    for q = 1:nl, fprintf(fid, ' %8.1f%%', 100*cnt(j,q)/max(s,1)); end
    fprintf(fid, '  (n=%d)\n', s);
end
fprintf(fid, '\n三个粗时段（扰动前 + 恢复晚期）合并：\n  ');
for q = 1:nl, fprintf(fid, '%s %.0f%%  ', labs{q}, 100*csum(q)/max(cn,1)); end
fprintf(fid, '(n=%d)\n', cn);
fclose(fid);
save(fullfile(outdir, 'event_labels.mat'), 'T', 'labs', 'cnt', 'edges', 'allab', 'atc');
fprintf('已写出：%s\n', fullfile(outdir, 'event_labels.txt'));
end
