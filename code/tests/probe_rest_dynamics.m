function probe_rest_dynamics()
% PROBE_REST_DYNAMICS  静息节点在不给任何外部触发时，增益能不能自己爬出去？
%
% 三区网络要"自发、再恢复、再等下一次"，前提是：静息点本身在慢动力学下不稳定，
% 即噪声驱动的 x = max(0, y_p - 0.03) 能把 G_e 一点点推上去。若不能，三区就永远
% 停在背景里，"自发"必须另找机制。
%
% 本脚本对比两种增益律（legacy / bistable），只给背景噪声，esf_enable = false，
% 跑 T 秒，输出：三个增益的首末值与极值、rms 的极大值、是否出现放电、
% 以及（若有）放电起点与下一次起点的间隔。
%
% 输出：tests/results/rest_dynamics.txt，并在 command window 打表。

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 600;
cases = { ...
    'legacy',   3.25, 22, 10; ...
    'legacy',   3.25, 35, 10; ...
    'bistable', 3.25, 22, 10; ...
    'bistable', 3.25, 35, 10};

fid = fopen(fullfile(resdir, 'rest_dynamics.txt'), 'w');
fprintf(fid, '# no ESF, single node, T=%d s; G values in mV\n', T);
fprintf(fid, '# law      A    B   G   Ge_end Gs_end Gf_end  Ge_max  Gs_min  Gf_max  rms_base rms_max  n_bursts  first_burst  period\n');

fprintf('\n%-9s %5s %4s %4s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n', ...
    'law','B','G','Ge_end','Gs_end','Gf_end','Ge_max','Gs_min','rms_base','rms_max','bursts','first');

for k = 1:size(cases, 1)
    law = cases{k,1};
    cfg = default_config();
    cfg.node_count = 1;
    cfg.A = cases{k,2}; cfg.B = cases{k,3}; cfg.G = cases{k,4};
    cfg.model_dynamic = 1;
    cfg.gain_law = law;
    cfg.esf_enable = false;      % 关键：不给任何外部触发
    cfg.simulation_time = T;
    cfg.resample_dy = true;
    cfg.seed = 20261101;

    [y, dy] = simulate_network(cfg);      % y: 1×n, dy: 3×n
    fs = cfg.eeg_fs;
    t = (0:numel(y)-1)/fs;

    Ge = dy(1, :); Gs = dy(2, :); Gf = dy(3, :);

    % rms 用 2 s 滑窗，判"放电"= 超过前 20 s 基线 rms 的 3 倍
    w = round(2*fs);
    nr = floor(numel(y)/w);
    rr = zeros(1, nr);
    for i = 1:nr
        seg = y((i-1)*w + (1:w));
        rr(i) = sqrt(mean(seg.^2));
    end
    tr = (0:nr-1)*2 + 1;
    base = median(rr(tr <= 20));
    thr  = 3*base;
    above = rr > thr;

    % 连段计数 = burst 数
    d = diff([0 above 0]);
    nb = sum(d == 1);
    on = tr(d(1:end-1) == 1);
    if numel(on) >= 2
        per = mean(diff(on));
    else
        per = NaN;
    end

    if isempty(on)
        tf = NaN;
    else
        tf = on(1);
    end

    str = sprintf('%-9s %5.0f %4.0f %8.2f %8.2f %8.2f %8.2f %8.2f %8.3f %8.3f %6d %8.1f %8.1f', ...
        law, cfg.B, cfg.G, Ge(end), Gs(end), Gf(end), max(Ge), min(Gs), base, max(rr), nb, tf, per);
    fprintf('%s\n', str);
    fprintf(fid, '%s\n', str);
end
fclose(fid);
fprintf('\n已写 results/rest_dynamics.txt\n');
end
