function probe_rest_x()
% PROBE_REST_X  静息点的"增益驱动"到底有多大：x = max(0, y_0 - 0.03)。
%
% 背景：probe_rest_dynamics 显示，从 (3.25,22,10) 出发、不给 ESF、跑 600 s，
% 三个增益一动不动（3.25→3.25，精确相等）。原因是驱动项取 max(0, y_0-0.03)，
% 而静息时锥体 PSP 的 y_0 低于 0.03，x 恒为 0，于是 dG 只剩恢复项 = 0。
%
% 三区网络若要"自发、恢复、等待下一次"，就必须让易感区的 x 偶尔为正。
% 本脚本扫 (A,B,G)，量 y_0 的均值、p95、以及 x>0 的时间占比（占空比）：
%   占空比 = 0    → 永远不会自发
%   占空比 1–5%   → 偶发，Ge 缓慢爬升（数十至数百秒），符合"等待下一次"
%   占空比 >20%   → 一上来就爬，等价于已经处在发作区
%
% 输出：tests/results/rest_x.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

Av = [3.25 3.60 4.00 4.60];
Bv = [10 15 20 25 30 45];
Gv = [5 10 25];
T  = 60;  T0 = 20;
dt = 1/256;  fs = 256;
y_thr = 0.03;

fid = fopen(fullfile(resdir, 'rest_x.txt'), 'w');
fprintf(fid, '# 静态增益，60 s（丢前 20 s）；x = max(0, y0-0.03)，y0 = 锥体 PSP 状态\n');
fprintf(fid, '#   A     B     G    y0_mean   y0_p95   x_mean   duty   rms_out   frac_over\n');

fprintf('\n%5s %5s %5s %9s %9s %9s %7s %9s\n', 'A','B','G','y0_mean','y0_p95','x_mean','duty','rms_out');
n = round(T/dt);
for ia = 1:numel(Av)
    for ib = 1:numel(Bv)
        for ig = 1:numel(Gv)
            cfg = default_config();
            cfg.A = Av(ia); cfg.B = Bv(ib); cfg.G = Gv(ig);
            m = wendling(cfg.A, cfg.B, cfg.G, dt, 0);   % 静态增益
            ns = gauss_noise(cfg, 1, n);
            y0v = zeros(1, n); outv = zeros(1, n);
            for i = 1:n
                [~, xi, yi] = m.simulate(ns(i));
                y0v(i) = xi(1);
                outv(i) = yi;
            end
            k = round(T0/dt)+1 : n;
            y0 = y0v(k); o = outv(k);
            x = max(0, y0 - y_thr);
            duty = mean(x > 0);
            fo   = mean(o > y_thr);
            ym = mean(y0); yp = prctile(y0, 95); xm = mean(x);
            rmsv = sqrt(mean(o.^2));

            str = sprintf('%5.2f %5.0f %5.0f %9.3f %9.3f %9.4f %7.3f %9.3f', ...
                Av(ia), Bv(ib), Gv(ig), ym, yp, xm, duty, rmsv);
            fprintf('%s\n', str);
            fprintf(fid, '%s\n', str);
        end
    end
end
fclose(fid);
fprintf('\n已写 results/rest_x.txt（共 %d 个参数点）\n', numel(Av)*numel(Bv)*numel(Gv));
end
