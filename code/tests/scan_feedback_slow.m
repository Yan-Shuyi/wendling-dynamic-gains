function scan_feedback_slow()
% SCAN_FEEDBACK_SLOW  把 PZ 调"迟钝"，给反向边一个可用的时间窗。
%
%   scan_feedback_slow()
%
% 为什么要这一组
%   scan_feedback 里 PZ→SOZ 权重加到 0.4 也没有结构性效果：SOZ 的放电幅度涨了 12%、
%   PZ 起振提前 3 s，但崩塌时刻、G_e 恢复、事件个数全不变。原因不是反向边没用，
%   而是**没有时间窗**：SOZ 与 PZ 的崩塌只差 2.7 s（188.7 vs 191.4），等 PZ 想回灌，
%   SOZ 自己也已经塌了。回路的延迟是 15 ms，事件的动力学尺度是几十秒，两者差三个数量级。
%
%   所以这里的对照组把 PZ 调迟钝（B 32→40、G 12→18），让它起振更晚、崩塌更晚，
%   在 SOZ 已经熄灭之后仍然在放电。如果"影响回来"是真的，这一组里 SOZ 的崩塌后
%   G_e 就不该回到 3.0，事件次数/间隔也该变。
%
% 输出
%   tests/results/feedback_slow.txt
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 900; C = 15;
SOZ = [3.25 22  8 0.005];
D   = [0 0 0; 0.015 0 0; 0.025 0 0];

PZsets = { ...
    'PZ_fast', [3.25 32 12 0]; ...   % 与 scan_feedback 相同（对照）
    'PZ_slow', [3.25 40 18 0]};      % 迟钝版
wlist = [0 0.40];
nm = {'SOZ','PZ','NIZ'};
NIZ = [3.25 50 22 0];

fid = fopen(fullfile(resdir, 'feedback_slow.txt'), 'w');
fprintf(fid, '# C=%d, T=%d s, 单段噪声；PZ 两组 × w_ret 两档\n', C, T);
fprintf(fid, '# pz        w_ret  zone  n_ev  t_ev1   t_ev2   Ge_min_after  rms_max\n');
fprintf('%-9s %-6s %-5s %5s %7s %7s %14s %8s\n', ...
    'pz','w_ret','zone','n_ev','t_ev1','t_ev2','Ge_col_end','rms_max');

for p = 1:size(PZsets,1)
    PZ = PZsets{p,2};
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
        for i = 1:3
            xi = highpass(y(i,:)', 0.5, cfg.eeg_fs)';
            segs = events_of(xi, fs);
            Ge = dy(3*(i-1)+1, :);
            if isempty(segs)
                gmin = NaN;
            else
                k1 = round((segs(1,1)+10)*fs) + 1;
                k2 = min(numel(Ge), round((segs(1,1)+60)*fs));
                gmin = min(Ge(k1:k2));
            end
            w = round(1*fs); nr = floor(numel(xi)/w);
            rr = zeros(1,nr);
            for j = 1:nr
                rr(j) = sqrt(mean(xi((j-1)*w+(1:w)).^2));
            end
            t1 = NaN; t2 = NaN;
            if size(segs,1) >= 1, t1 = segs(1,1); end
            if size(segs,1) >= 2, t2 = segs(2,1); end
            fprintf('%-9s %-6.2f %-5s %5d %7.1f %7.1f %14.3f %8.2f\n', ...
                PZsets{p,1}, wlist(wi), nm{i}, size(segs,1), t1, t2, gmin, max(rr));
            fprintf(fid, '%-9s %-6.2f %-5s %5d %7.1f %7.1f %14.3f %8.2f\n', ...
                PZsets{p,1}, wlist(wi), nm{i}, size(segs,1), t1, t2, gmin, max(rr));
        end
    end
end
fclose(fid);
fprintf('\n已写 results/feedback_slow.txt\n');
end

function segs = events_of(x, fs)
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
        j = k + round(20/0.25);
    else
        j = j + 1;
    end
end
end
