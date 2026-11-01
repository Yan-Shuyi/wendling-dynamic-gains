function analyze_three_zone()
% ANALYZE_THREE_ZONE  从 three_zone.mat 重算跨区次序（不重跑仿真）。
%
%   analyze_three_zone()
%
% 为什么单独一步
%   run_three_zone 里的 onset 用的是"3×基线 / 1 s 窗"，而在三区里基线差了 3 倍
%   （SOZ 1.00、PZ 0.31、NIZ 0.38），同一个相对阈值在不同区等价于不同的绝对门槛，
%   于是 PZ 只要被耦合轻推一下就"先起"。这里把判据拆成两层：
%     t_elec  电学跟随：0.25 s 窗 rms 首次持续 2 s 越过 max(4×基线, 0.5 mV)
%             （0.5 mV 只是防死平的兜底，不是幅度门槛；去直流后三区基线一致，
%             4× 基线就等于同一个绝对门槛）
%     t_gain  增益募集：G_e 首次超过 3.60 mV（静息 3.25，驱动阈值约 3.8）
%   前者看"同一事件有没有立刻传过去"，后者看"对方自己有没有被点着"。
%
% 输出
%   tests/results/three_zone_onset.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');
L = load(fullfile(resdir, 'three_zone.mat'));

fs = 256;
names = {'SOZ','PZ','NIZ'};
tags  = fieldnames(L);

fid = fopen(fullfile(resdir, 'three_zone_onset.txt'), 'w');
fprintf(fid, '# 跨区次序（两层判据）：t_elec = 电学跟随，t_gain = 增益募集\n');
fprintf(fid, '# case        zone   rms_base  thr    t_elec   t_gain3.6  t_gain4.0  Ge_max  rms_burst\n');

for k = 1:numel(tags)
    D = L.(tags{k});
    n = size(D.y, 1);
    fprintf('\n=== %s ===\n', tags{k});
    fprintf('%-6s %8s %6s %9s %9s %9s %8s %9s\n', ...
        'zone','rms_b','thr','t_elec','t_Ge>3.6','t_Ge>4.0','Ge_max','rms_burst');

    for i = 1:n
        % 0.25 s 窗 rms（步长 0.25 s），先去 0.5 Hz 高通
        xi = highpass(D.y(i,:)', 0.5, fs)';
        w = round(0.25*fs); nr = floor(numel(xi)/w);
        rr = zeros(1, nr);
        for j = 1:nr
            seg = xi((j-1)*w+(1:w));
            rr(j) = sqrt(mean(seg.^2));
        end
        tr = (0:nr-1)*0.25 + 0.125;
        base = median(rr(tr <= 30));
        thr  = max(4*base, 0.5);
        above = rr > thr;
        need = 8;                     % 持续 2 s（8 × 0.25 s）
        te = NaN;
        for j = 1:nr-need
            if all(above(j:j+need-1)), te = tr(j) - 0.125; break; end
        end

        Ge = D.dy(3*(i-1)+1, :); t = D.t;
        g36 = first_cross(Ge, t, 3.60);
        g40 = first_cross(Ge, t, 4.00);
        rms_burst = max(rr(tr > 150 & tr < 700));

        fprintf('%-6s %8.3f %6.1f %9.1f %9.1f %9.1f %8.2f %9.2f\n', ...
            names{i}, base, thr, te, g36, g40, max(Ge), rms_burst);
        fprintf(fid, '%-12s %-5s %8.3f %6.1f %9.1f %9.1f %9.1f %8.2f %9.2f\n', ...
            tags{k}, names{i}, base, thr, te, g36, g40, max(Ge), rms_burst);
    end
end
fclose(fid);

%% ---- 电学传播的毫秒级滞后（SOZ 放电窗内）----
D = L.net;
t0 = D.onset(1) + 30;                       % 放电已稳定
k1 = round(t0*fs) + 1; k2 = min(size(D.y,2), k1 + round(120*fs) - 1);
x0 = D.y(1, k1:k2); x0 = x0 - mean(x0);
maxlag = round(0.05*fs);                    % ±50 ms
fprintf('\n--- 与 SOZ 的互相关峰（%.0f–%.0f s，±50 ms）---\n', t0, t0+120);
for i = 2:3
    y0 = D.y(i, k1:k2); y0 = y0 - mean(y0);
    [c, lg] = xcorr(y0, x0, maxlag, 'coeff');
    [cm, im] = max(c);
    fprintf('%-4s  lag = %+6.1f ms   r = %.3f\n', names{i}, lg(im)/fs*1000, cm);
end

%% ---- 崩塌时刻：Ge 峰值与 Gs 谷值 ----
n = 3;
fprintf('\n--- 第一周期刻度（s）---\n');
fprintf('%-6s %10s %10s %10s %10s\n', 'zone','Ge_peak','Gs_trough','Gs<10','Ge>6');
for i = 1:n
    Ge = D.dy(3*(i-1)+1, :); Gs = D.dy(3*(i-1)+2, :); Gf = D.dy(3*(i-1)+3, :);
    [~, ip] = max(Ge);
    [~, it] = min(Gs);
    fprintf('%-6s %10.1f %10.1f %10.1f %10.1f\n', names{i}, ...
        D.t(ip), D.t(it), first_cross(-Gs, D.t, -10), first_cross(Ge, D.t, 6.0));
end
end

function t0 = first_cross(x, t, thr)
k = find(x > thr, 1);
if isempty(k), t0 = NaN; else, t0 = t(k); end
end
