function probe_duty_fine()
% PROBE_DUTY_FINE  在 B = 15–20 mV 之间细扫：静息点"偶发"的窗口有多宽。
%
% probe_rest_x 显示占空比（x>0 的时间比例）在 B 上跳变：B=15 时 duty=1，
% B=20 时 duty=0，中间没有过渡点。而三区若要"缓慢自发"，需要 1%–10% 的占空比。
% 本脚本把 B 步长缩到 0.25 mV，A 取 3.25/3.60，G 取 10，找这个窗口是否存在、有多宽。
%
% 输出：tests/results/duty_fine.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

Bv = 15:0.25:20;
Av = [3.25 3.60];
G  = 10;
T = 60; T0 = 20; dt = 1/256; y_thr = 0.03;

fid = fopen(fullfile(resdir, 'duty_fine.txt'), 'w');
fprintf(fid, '# A=3.25/3.60, G=10, static gains, 60 s (drop 20 s)\n');
fprintf(fid, '#   A      B   y0_mean   y0_p95   x_mean    duty   rms_out\n');
for ia = 1:numel(Av)
    for ib = 1:numel(Bv)
        cfg = default_config();
        cfg.A = Av(ia); cfg.B = Bv(ib); cfg.G = G;
        m = wendling(cfg.A, cfg.B, cfg.G, dt, 0);
        ns = gauss_noise(cfg, 1, round(T/dt));
        n = numel(ns);
        y0v = zeros(1,n); outv = zeros(1,n);
        for i = 1:n
            [~, xi, yi] = m.simulate(ns(i));
            y0v(i) = xi(1); outv(i) = yi;
        end
        k = round(T0/dt)+1 : n;
        y0 = y0v(k); o = outv(k);
        x = max(0, y0 - y_thr);
        str = sprintf('%6.2f %6.2f %8.3f %8.3f %8.4f %7.4f %8.3f', ...
            Av(ia), Bv(ib), mean(y0), prctile(y0,95), mean(x), mean(x>0), sqrt(mean(o.^2)));
        fprintf('%s\n', str); fprintf(fid, '%s\n', str);
    end
end
fclose(fid);
fprintf('\n已写 results/duty_fine.txt\n');
end
