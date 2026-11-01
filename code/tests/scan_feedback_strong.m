function scan_feedback_strong()
% SCAN_FEEDBACK_STRONG  反向边加到与前馈同量级（0.6–1.0）：看"影响回来"的边界。
%
%   scan_feedback()
%
% 背景
%   run_three_zone 的 net 是纯前馈级联（只有 SOZ→PZ、SOZ→NIZ），所以 PZ 被募集
%   之后对 SOZ 没有任何作用，NIZ 更是只读。要谈"互相影响"，最小改动是加一条
%   PZ→SOZ 的边。这里扫 w = 0 / 0.1 / 0.2 / 0.4，其余全部不动（C=15、单段噪声、同种子）。
%
% 看什么
%   1. SOZ 的崩塌时刻是否提前、G_s 谷值是否更深（被 PZ 的输入压下去）；
%   2. SOZ 崩塌后 G_e 是否回不到 3.25（被 PZ 续着），从而不需要等 bias 慢慢爬；
%   3. 第二个事件的起点是否提前 → 周期是否缩短（这是"影响回来"最硬的指标：
%      网络的回程改变了 SOZ 自己的节律，而不只是加了一点噪声）。
%
% 输出
%   tests/results/feedback.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 900; C = 15;
SOZ = [3.25 22  8 0.005];
PZ  = [3.25 32 12 0];
NIZ = [3.25 50 22 0];
D   = [0 0 0; 0.015 0 0; 0.025 0 0];

wlist = [0 0.60 0.80 1.00];
nm = {'SOZ','PZ','NIZ'};

fid = fopen(fullfile(resdir, 'feedback_strong.txt'), 'w');
fprintf(fid, '# C=%d, T=%d s, 单段噪声, 只改 W(SOZ<-PZ)\n', C, T);
fprintf(fid, '# w_ret  zone  n_ev  t_ev1   t_ev2   Ge_min_after  Gs_nadir  rms_max\n');
fprintf('%-6s %-5s %5s %7s %7s %14s %9s %8s\n', ...
    'w_ret','zone','n_ev','t_ev1','t_ev2','Ge_collapse_end','Gs_nadir','rms_max');

S = struct();
for wi = 1:numel(wlist)
    W = [0 0 0; 1 0 0; 0.15 0 0];
    if wlist(wi) > 0, W(1,2) = wlist(wi); end
    cfg = default_config();
    cfg.node_count = 3;
    cfg.A = [SOZ(1) PZ(1) NIZ(1)];
    cfg.B = [SOZ(2) PZ(2) NIZ(2)];
    cfg.G = [SOZ(3) PZ(3) NIZ(3)];
    cfg.gain_bias = [SOZ(4) PZ(4) NIZ(4)];
    cfg.model_dynamic = 1; cfg.gain_law = 'legacy'; cfg.gain_reset_tol = 0.05;
    cfg.esf_enable = false; cfg.simulation_time = T; cfg.resample_dy = true;
    cfg.seed = 20261101; cfg.coupling_gain = C;
    cfg.sc_matrix = W; cfg.delay_matrix = D;
    cfg.noise_segments = 1;

    [y, dy, cfg] = simulate_network(cfg);
    fs = cfg.eeg_fs;

    ev = cell(1,3);
    for i = 1:3
        xi = highpass(y(i,:)', 0.5, fs)';   % 事件判据必须在去直流后的信号上，与 rms 同口径
        [te, segs] = events_of(xi, fs);
        ev{i} = segs;
        Ge = dy(3*(i-1)+1, :);
        Gs = dy(3*(i-1)+2, :);
        % 第一次崩塌之后 10–60 s 的 G_e 最低值：若回不到静息 3.25，说明被对方续着
        if isempty(segs)
            gmin = NaN; gsmin = NaN;
        else
            k1 = round((segs(1,1) + 10)*fs) + 1;
            k2 = min(numel(Ge), round((segs(1,1) + 60)*fs));
            gmin = min(Ge(k1:k2));
            gsmin = min(Gs(1:min(numel(Gs), round((segs(1,1)+120)*fs))));
        end
        w = round(1*fs); nr = floor(numel(xi)/w);
        rr = zeros(1,nr);
        for j = 1:nr
            seg = xi((j-1)*w+(1:w));
            rr(j) = sqrt(mean(seg.^2));
        end
        t1 = NaN; t2 = NaN;
        if size(segs,1) >= 1, t1 = segs(1,1); end
        if size(segs,1) >= 2, t2 = segs(2,1); end
        fprintf('%-6.2f %-5s %5d %7.1f %7.1f %14.3f %9.2f %8.2f\n', ...
            wlist(wi), nm{i}, size(segs,1), t1, t2, gmin, gsmin, max(rr));
        fprintf(fid, '%-6.2f %-5s %5d %7.1f %7.1f %14.3f %9.2f %8.2f\n', ...
            wlist(wi), nm{i}, size(segs,1), t1, t2, gmin, gsmin, max(rr));
    end
    e.y = y; e.dy = dy; e.t = (0:size(y,2)-1)/fs; e.ev = ev; e.w = wlist(wi);
    S.(sprintf('w%03d', round(wlist(wi)*100))) = e;
end
fclose(fid);
save(fullfile(resdir, 'feedback_strong.mat'), '-struct', 'S');
fprintf('\n已写 results/feedback.txt / .mat\n');
end

% ----------------------------------------------------------------------
function [t0, segs] = events_of(x, fs)
% 0.25 s 窗 rms，持续 2 s 越过 max(4×基线, 2.0 mV) 记为一个事件；段间至少隔 20 s
w = round(0.25*fs); nr = floor(numel(x)/w);
rr = zeros(1,nr);
for j = 1:nr
    rr(j) = sqrt(mean(x((j-1)*w+(1:w)).^2));
end
tr = (0:nr-1)*0.25 + 0.125;
base = median(rr(tr <= 30));
above = rr > max(4*base, 0.5);
segs = zeros(0,2);
j = 1;
while j <= nr-8
    if all(above(j:j+7))
        k = j + 7;
        while k < nr && any(above(k:min(nr,k+3))), k = k + 1; end
        segs(end+1,:) = [tr(j)-0.125, tr(k)-0.125]; %#ok<AGROW>
        j = k + round(20/0.25);          % 跳过 20 s，避免一个事件被切碎
    else
        j = j + 1;
    end
end
t0 = tr; %#ok<NASGU>
end
