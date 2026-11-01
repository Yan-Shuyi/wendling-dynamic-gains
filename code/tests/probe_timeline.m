function probe_timeline()
% PROBE_TIMELINE  逐 5 s 窗列出 fig2 配置的幅度/主频/频带占比与三个增益，
% 用来判定"EEG 上真正的事件时刻"和"增益机器的事件时刻"是否一致。
%
% 输出：results/timeline.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();
p2 = classify_params_v2();

cfg = default_config();
cfg.node_count = 1;
cfg.esf_enable = true; cfg.esf_time = 10; cfg.esf_params = [4.25 22 10];
[out, dy, cfg] = simulate_network(cfg);
fs = cfg.eeg_fs; dfs = cfg.dy_fs;

ref = make_reference(out(1:min(numel(out), 9*fs)), fs, p2);

fid = fopen(fullfile(outdir, 'timeline.txt'), 'w');
fprintf(fid, 'fig2 配置：esf at 10 s, esf_params=[4.25 22 10]\n');
fprintf(fid, '%-7s %7s %7s %6s %6s %6s %7s %7s %7s  %s\n', ...
    't(win)', 'rms', 'f_peak', 'slow%', '4-10%', '15-40%', 'Ge', 'Gs', 'Gf', 'v2 判定');
for t0 = 0:5:195
    i1 = max(1, round(t0*fs)+1); i2 = min(numel(out), round((t0+5)*fs));
    x = out(i1:i2);
    f = extract_features_v2(x, fs, p2);
    [lab, sub, ~, sc] = classify_v2(x, fs, p2, ref);
    k = min(numel(dy), round((t0+2.5)*dfs)+1);
    if isempty(sub), s = lab; else, s = [lab '/' sub]; end
    fprintf(fid, '%3.0f-%3.0f %7.2f %7.2f %6.1f %6.1f %6.1f %7.2f %7.2f %7.2f  %-16s rms_ratio=%.1f spk=%.2fHz\n', ...
        t0, t0+5, f.rms, f.f_peak, 100*f.slow, 100*f.rel_slow, 100*f.rel_fast, ...
        dy(1,k), dy(2,k), dy(3,k), s, sc.rms_ratio, sc.spk_rate);
end
fclose(fid);
end
