function probe_events()
% PROBE_EVENTS  把"EEG 上真的发生了什么"和"增益机器的事件"分开量：
%   t_EEG  = 1 s 滑窗 rms 首次超过 3×基线（电发作起点）
%   t_Ge   = 高幅段结束（rms 回落到 1.5×基线以下）
%   t_Gs   = Gs 首次跌破 12.75 mV（旧 D14 用的锚点）
%
% 输出：results/events.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();
fid = fopen(fullfile(outdir, 'events.txt'), 'w');
fprintf(fid, '触发固定 at 10 s；只改 esf 的 Ge\n');
fprintf(fid, '%6s %9s %9s %9s %9s %9s\n', 'Ge', 't_EEG', 't_high_end', 't_Gs<12.75', 'dur_EEG', 'dur_Gs');
for ge = [4.25 4.10 4.00 3.90]
    cfg = default_config(); cfg.node_count = 1;
    cfg.esf_enable = true; cfg.esf_time = 10; cfg.esf_params = [ge 22 10];
    [out, dy, cfg] = simulate_network(cfg);
    fs = cfg.eeg_fs; dfs = cfg.dy_fs;
    base = std(out(1:min(numel(out), 9*fs)));
    w = fs;                                  % 1 s 滑窗
    n = floor(numel(out)/w);
    r = zeros(1,n);
    for k = 1:n
        r(k) = std(out((k-1)*w+1 : k*w));
    end
    tk = (1:n);
    k1 = find(r > 3*base, 1, 'first');
    if isempty(k1), t_eeg = NaN; else, t_eeg = tk(k1)-1; end
    if isnan(t_eeg)
        t_end = NaN;
    else
        k2 = find(r(k1:end) < 1.5*base, 1, 'first');
        if isempty(k2), t_end = NaN; else, t_end = tk(k1+k2-1)-1; end
    end
    gs = dy(2,:);
    kg = find(gs < 12.75, 1, 'first');
    if isempty(kg), t_gs = NaN; else, t_gs = (kg-1)/dfs; end
    fprintf(fid, '%6.2f %9.1f %9.1f %9.1f %9.1f %9.1f\n', ge, t_eeg, t_end, t_gs, ...
        t_end - t_eeg, t_gs - 10);
end
fclose(fid);
end
