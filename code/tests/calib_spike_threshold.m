function calib_spike_threshold()
% CALIB_SPIKE_THRESHOLD  校准棘波斜率阈值：参考段用"白噪声"和"1/f 样 EEG"，看低通后
% 斜率分布如何、各原型棘波串的最大斜率落在哪，据此定 p.spk_slope_k 与 p.spk_lp。
%
% 同时量一下模型 ICTAL 段的主频，判断是否需要 delta 亚类。
%
% 输出：results/calib_spike.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();
fid = fopen(fullfile(outdir, 'calib_spike.txt'), 'w');
Fs = 256; T = 10; t = (0:1/Fs:T-1/Fs)';

rng(11);
% 三种参考段
refs = {};
refs(end+1,:) = {'白噪声 sigma=0.8', 0.8*randn(size(t))};
refs(end+1,:) = {'白噪声 sigma=2.0', 2.0*randn(size(t))};
pk = cumsum(0.5*randn(size(t)));  pk = pk - mean(pk);  pk = 1.2*pk/std(pk);
refs(end+1,:) = {'1/f 样 sigma=1.2', pk};

% 原型棘波串
spk = {};
for rt = [0.3 1 2 4]
    for amp = [4 8]
        spk(end+1,:) = {sprintf('棘波 %.1f Hz amp=%g', rt, amp), spike_train(t, rt, 0.02, amp)};
    end
end

lps = [0 40 30];
for il = 1:numel(lps)
    lp = lps(il);
    fprintf(fid, '\n========== 低通截止 lp = %d Hz ==========\n', lp);
    fprintf(fid, '参考段斜率分布（mV/s，低通后逐点差分×Fs）\n');
    fprintf(fid, '%-24s %8s %8s %8s\n', '参考段', 'p50', 'p99', 'max');
    base = zeros(1, numel(refs(:,1)));
    for r = 1:size(refs,1)
        d = slope_of(refs{r,2}, Fs, lp);
        base(r) = prctile(d, 99);
        fprintf(fid, '%-24s %8.1f %8.1f %8.1f\n', refs{r,1}, prctile(d,50), prctile(d,99), max(d));
    end
    fprintf(fid, '\n原型棘波最大斜率（mV/s）与"相对白噪声 p99 的倍数"\n');
    fprintf(fid, '%-24s %8s %10s %10s\n', '原型', 'max', '/白0.8p99', '/1/f p99');
    for s = 1:size(spk,1)
        d = slope_of(spk{s,2}, Fs, lp);
        fprintf(fid, '%-24s %8.1f %10.2f %10.2f\n', spk{s,1}, max(d), max(d)/base(1), max(d)/base(3));
    end
end

%% 模型 ICTAL 段主频（判断是否需要 delta 亚类）
fprintf(fid, '\n========== 模型分期的实测主频 ==========\n');
p2 = classify_params_v2();
cfg = default_config(); cfg.node_count = 1;
cfg.esf_enable = true; cfg.esf_time = 10; cfg.esf_params = [4.25 22 10];
[out, ~, cfg] = simulate_network(cfg);
fs = cfg.eeg_fs;
bins = {'interICTAL',5,10; 'preONSET',30,40; 'ONSET',46,49; 'ICTAL',55,65; 'post',80,120; 'late',150,190};
fprintf(fid, '%-12s %7s %7s %7s %7s %7s\n', '窗', 'f_peak', 'slow%', '4-10%', '15-40%', 'rms');
for b = 1:size(bins,1)
    i1 = max(1, round(bins{b,2}*fs)+1); i2 = min(numel(out), round(bins{b,3}*fs));
    f = extract_features_v2(out(i1:i2), fs, p2);
    fprintf(fid, '%-12s %7.2f %7.1f %7.1f %7.1f %7.2f\n', bins{b,1}, f.f_peak, ...
        100*f.slow, 100*f.rel_slow, 100*f.rel_fast, f.rms);
end
fclose(fid);
end

function d = slope_of(x, Fs, lp)
x = double(x(:)); x = x - mean(x);
if lp > 0
    [b,a] = butter(4, lp/(Fs/2));
    x = filtfilt(b, a, x);
end
d = abs(diff(x)) * Fs;
end

function x = spike_train(t, rate, width, amp)
x = zeros(size(t));
for c0 = 0.1:(1/rate):(t(end)-0.1)
    x = x + amp*exp(-((t-c0)/width).^2);
end
end
